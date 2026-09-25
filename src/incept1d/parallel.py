# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Running independent solves in parallel.

The points of an inception-curve sweep are independent: each is a root
search at its own p·d.  :func:`parallel_map` runs such tasks in worker
processes created by ``fork``, so the task — typically a closure over a
loaded mechanism, which is an exec'd module and does not pickle reliably —
is inherited by the workers rather than sent to them.  Only the task index
goes out and only the (picklable) result comes back.  Where ``fork`` is not
available the tasks run in this process instead.

:func:`physical_cores` is the default worker count: one worker per physical
core, never one per hardware thread, because two workers sharing a core
through hyperthreading only slow each other down on this workload.
"""

import multiprocessing
import os
import queue
import subprocess
import traceback


def _worker(task, blocks, results):
    """
    Body of one worker process: solve its blocks, report every point.

    Each block is solved in order, every point seeded by the one before.
    Results go back as ``(i, result)``; a failure as ``(None, traceback)``,
    so the parent can raise it instead of waiting for results that will
    never come.
    """
    try:
        for block in blocks:
            previous = None
            for i in block:
                previous = task(i, previous)
                results.put((i, previous))
    except BaseException:
        results.put((None, traceback.format_exc()))


def physical_cores():
    """
    The number of physical CPU cores, or 1 if it cannot be determined.

    Hardware threads are not counted: an 8+8 core machine with
    hyperthreading on its 8 performance cores has 24 logical CPUs but 16
    physical cores.  Read from psutil if installed, else from
    ``/proc/cpuinfo`` (Linux) or ``sysctl hw.physicalcpu`` (macOS); any
    other platform, or any failure, gives 1 rather than a guess that might
    put two workers on one core.
    """
    try:
        import psutil

        n = psutil.cpu_count(logical=False)
        if n:
            return int(n)
    except ImportError:
        pass
    try:
        cores, phys, core = set(), None, None
        with open("/proc/cpuinfo") as fh:
            for line in fh:
                key, _, value = line.partition(":")
                key = key.strip()
                if key == "physical id":
                    phys = value.strip()
                elif key == "core id":
                    core = value.strip()
                elif not key:
                    if phys is not None and core is not None:
                        cores.add((phys, core))
                    phys = core = None
        if phys is not None and core is not None:
            cores.add((phys, core))
        if cores:
            return len(cores)
    except OSError:
        pass
    try:
        out = subprocess.run(
            ["sysctl", "-n", "hw.physicalcpu"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        n = int(out.stdout.strip())
        if n > 0:
            return n
    except (OSError, ValueError, subprocess.SubprocessError):
        pass
    return 1


def parallel_map(task, n, jobs, on_result=None, blocks_per_job=1):
    """
    Evaluate ``task(i, previous)`` for ``i`` in ``range(n)``, in parallel.

    The indices are split into contiguous blocks, each solved in order by
    one worker, and *previous* is the result of the index before it in the
    same block (None for the first).  That keeps what makes a sweep cheap —
    each point warm-started from its neighbour — while the blocks run
    concurrently; only the first point of each block starts cold.  A
    shared queue of single indices would give every point a cold start
    instead, which on a strongly non-uniform gap costs more than twice the
    CPU time.

    Parameters
    ----------
    task : callable
        ``task(i, previous) -> result``.  Need not be picklable; the result
        must be.
    n : int
        Number of tasks.
    jobs : int
        Maximum number of worker processes.  1, a single task, or a
        platform without ``fork`` runs everything in this process, in one
        block.
    on_result : callable or None
        Called as ``on_result(i, result)`` in this process as each task
        finishes — in completion order, which differs from index order when
        running in parallel.
    blocks_per_job : int
        Blocks per worker.  More, smaller blocks balance uneven costs better
        at the price of one cold start each.

    Returns
    -------
    list
        The results, in index order.
    """
    results = [None] * n
    workers = min(int(jobs), n)
    if workers <= 1 or "fork" not in multiprocessing.get_all_start_methods():
        previous = None
        for i in range(n):
            previous = results[i] = task(i, previous)
            if on_result is not None:
                on_result(i, previous)
        return results

    n_blocks = min(n, workers * max(1, int(blocks_per_job)))
    edges = [round(k * n / n_blocks) for k in range(n_blocks + 1)]
    blocks = [range(edges[k], edges[k + 1]) for k in range(n_blocks)]

    # Plain forked processes, all started before this process creates any
    # thread of its own.  A pool that forks on demand does so while its
    # management thread runs, and a lock held by any thread at that moment
    # stays held forever in the child: under pytest that deadlocked.
    context = multiprocessing.get_context("fork")
    results_q = context.Queue()
    procs = [
        context.Process(
            target=_worker, args=(task, blocks[w::workers], results_q), daemon=True
        )
        for w in range(workers)
    ]
    for proc in procs:
        proc.start()
    try:
        received = 0
        while received < n:
            try:
                i, result = results_q.get(timeout=1.0)
            except queue.Empty:
                if not any(proc.is_alive() for proc in procs):
                    raise RuntimeError(
                        f"worker processes exited after {received} of {n} results"
                    )
                continue
            if i is None:
                raise RuntimeError(f"a worker process failed:\n{result}")
            results[i] = result
            received += 1
            if on_result is not None:
                on_result(i, result)
    finally:
        for proc in procs:
            if proc.is_alive():
                proc.terminate()
            proc.join()
    return results


def limit_blas_threads():
    """
    Ask the BLAS and OpenMP runtimes for one thread each.

    Must run before NumPy is imported to take effect.  The matrices here
    are tiny (at most a few dozen rows), and a multithreaded BLAS spreads
    them over many cores for no gain: a pdiv sweep measured 16.9 s at
    ~1200 % CPU with the default and 14.6 s at 100 % with one thread.  It
    also leaves the cores to the worker processes of :func:`parallel_map`.
    An explicit setting in the environment is left alone.
    """
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "BLIS_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ.setdefault(name, "1")
