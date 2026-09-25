# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Parallel sweeps: same answers as sequential ones, in the same order."""

import functools
import os
import time

import numpy as np
import pytest

from incept1d.inception import compute_inception_curve
from incept1d.parallel import parallel_map, physical_cores
from incept1d.solver import inception_det


def test_physical_cores_is_a_positive_count():
    n = physical_cores()
    assert isinstance(n, int) and n >= 1


class TestParallelMap:
    def test_results_come_back_in_index_order(self):
        assert parallel_map(lambda i, prev: i * i, 20, jobs=4) == [
            i * i for i in range(20)
        ]

    def test_task_need_not_be_picklable(self):
        """Workers are forked, so a closure over local state is fine."""
        offset = {"value": 7}
        assert parallel_map(lambda i, prev: i + offset["value"], 5, jobs=3) == [
            7,
            8,
            9,
            10,
            11,
        ]

    def test_every_result_is_reported_once(self):
        seen = []
        parallel_map(
            lambda i, prev: -i, 12, jobs=4, on_result=lambda i, r: seen.append((i, r))
        )
        assert sorted(seen) == [(i, -i) for i in range(12)]

    def test_blocks_are_contiguous_and_chained(self):
        """Each index sees its predecessor's result, except at block starts."""
        out = parallel_map(
            lambda i, prev: (i, None if prev is None else prev[0]), 12, jobs=3
        )
        starts = [i for i, (_, prev) in enumerate(out) if prev is None]
        assert starts == [0, 4, 8]
        assert all(prev == i - 1 for i, (_, prev) in enumerate(out) if prev is not None)

    def test_one_dead_worker_is_reported_while_others_run(self):
        """
        A worker that dies without a traceback (killed, or a crash in native
        code) owes results it will never send.  The sweep must fail promptly,
        not wait for them while another worker is still busy.
        """

        def task(i, prev):
            if i == 0:
                os._exit(1)
            time.sleep(60.0)
            return i

        start = time.monotonic()
        with pytest.raises(RuntimeError, match="exited after 0 of 2"):
            parallel_map(task, 2, jobs=2)
        assert time.monotonic() - start < 15.0

    def test_one_job_runs_in_this_process(self):
        order = []
        parallel_map(lambda i, prev: order.append(i), 5, jobs=1)
        assert order == [0, 1, 2, 3, 4]


class TestParallelCurve:
    """jobs > 1 must give the same branches as jobs = 1."""

    @pytest.mark.parametrize("all_branches", [False, True])
    def test_same_branches(self, toy, uniform, all_branches):
        pd_arr = np.geomspace(1e-3, 1.0, 12)
        det = functools.partial(inception_det, field_dist=uniform)
        runs = [
            compute_inception_curve(
                pd_arr, toy, 1.0, 293.0, all_branches=all_branches, det_fn=det, jobs=j
            )
            for j in (1, 3)
        ]
        assert len(runs[0]) == len(runs[1]) > 0
        for a, b in zip(*runs):
            assert np.array_equal(a["idx"], b["idx"])
            assert a["EN"] == pytest.approx(b["EN"], rel=1e-9)

    def test_progress_reports_every_point(self, toy, uniform):
        pd_arr = np.geomspace(1e-3, 1.0, 8)
        seen = []
        compute_inception_curve(
            pd_arr,
            toy,
            1.0,
            293.0,
            det_fn=functools.partial(inception_det, field_dist=uniform),
            progress=lambda i, n, pd, p, roots, t, check: seen.append((i, n)),
            jobs=3,
        )
        assert sorted(seen) == [(i, 8) for i in range(8)]
