# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Safe arithmetic expressions for declarative mechanism files.

A mechanism file writes rate coefficients, mobilities and yields as formulas,
``"1.24e-17*exp(-(259/(14.6 + EN))**2)"``.  This module turns such a string
into a callable without giving the file the power of Python: the string is
parsed with :mod:`ast`, every node of the tree is checked against a short
whitelist — numbers, names, the arithmetic operators and a few mathematical
functions — and only a tree that passes is compiled.  There are no
attribute accesses, subscripts, comparisons, conditionals, lambdas,
comprehensions or calls to anything but the listed functions, so an
expression can compute a number and nothing else.

Expressions are evaluated with NumPy scalars, so a division by zero gives
``inf`` (and ``exp(-inf) = 0``) rather than an exception, as the Python
mechanisms relied on.
"""

import ast

import numpy as np

#: Functions an expression may call, by name.
FUNCTIONS = {
    "exp": np.exp,
    "log": np.log,
    "log10": np.log10,
    "sqrt": np.sqrt,
    "abs": np.abs,
    "tanh": np.tanh,
    "cos": np.cos,
    "sin": np.sin,
    "min": np.minimum,
    "max": np.maximum,
}

_BINARY = (ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow)
_UNARY = (ast.UAdd, ast.USub)


class ExpressionError(ValueError):
    """An expression that does not parse, uses what is not allowed, or names
    an unknown variable."""


class Expression:
    """
    A checked, compiled formula.

    Parameters
    ----------
    text : str or float
        The formula.  A number is accepted as a constant.
    names : iterable of str
        The variables the formula may use.
    where : str
        Where the formula comes from, for error messages.

    Attributes
    ----------
    text : str
        The formula as written.
    names : frozenset of str
        The variables it actually uses.
    """

    def __init__(self, text, names, where="expression"):
        if isinstance(text, bool) or not isinstance(text, (str, int, float)):
            raise ExpressionError(
                f"{where}: expected a formula or a number, got {text!r}"
            )
        self.text = str(text)
        try:
            tree = ast.parse(self.text.strip(), mode="eval")
        except SyntaxError as exc:
            raise ExpressionError(
                f"{where}: cannot parse {self.text!r}: {exc.msg}"
            ) from None
        allowed = set(names)
        used = set()
        callees = {id(n.func) for n in ast.walk(tree) if isinstance(n, ast.Call)}
        for node in ast.walk(tree):
            _check(node, allowed, used, callees, where, self.text)
        self.names = frozenset(used)
        # Every literal becomes a NumPy float bound by name: Python integer
        # arithmetic such as 9**9**9 is exact and unbounded, Python float
        # arithmetic raises on overflow, NumPy's overflows to inf.
        self._literals = {}
        tree = _Literals(self._literals).visit(tree)
        self._code = compile(ast.fix_missing_locations(tree), f"<{where}>", "eval")

    def __call__(self, env):
        """Evaluate with the variables in *env* (a mapping name -> value)."""
        scope = {k: np.float64(env[k]) for k in self.names}
        scope.update(self._literals)
        scope.update(FUNCTIONS)
        with np.errstate(divide="ignore", over="ignore"):
            return float(eval(self._code, {"__builtins__": {}}, scope))  # noqa: S307

    def __repr__(self):
        return f"Expression({self.text!r})"


class _Literals(ast.NodeTransformer):
    """Replace each numeric literal by a name bound to its NumPy float."""

    def __init__(self, literals):
        self.literals = literals

    def visit_Constant(self, node):
        name = f"_c{len(self.literals)}"  # not a valid variable: names are checked
        self.literals[name] = np.float64(node.value)
        return ast.copy_location(ast.Name(id=name, ctx=ast.Load()), node)


def _check(node, allowed, used, callees, where, text):
    """Raise ExpressionError unless *node* is on the whitelist."""
    if isinstance(node, (ast.Expression, ast.Load)):
        return
    if isinstance(node, ast.Constant):
        if isinstance(node.value, bool) or not isinstance(node.value, (int, float)):
            raise ExpressionError(f"{where}: only numbers are allowed in {text!r}")
        return
    if isinstance(node, ast.BinOp):
        if not isinstance(node.op, _BINARY):
            raise ExpressionError(
                f"{where}: operator {type(node.op).__name__} is not allowed in {text!r}"
            )
        return
    if isinstance(node, ast.UnaryOp):
        if not isinstance(node.op, _UNARY):
            raise ExpressionError(
                f"{where}: operator {type(node.op).__name__} is not allowed in {text!r}"
            )
        return
    if isinstance(node, (ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow)):
        return
    if isinstance(node, (ast.UAdd, ast.USub)):
        return
    if isinstance(node, ast.Call):
        if (
            not isinstance(node.func, ast.Name)
            or node.func.id not in FUNCTIONS
            or node.keywords
        ):
            raise ExpressionError(
                f"{where}: only {sorted(FUNCTIONS)} may be called, in {text!r}"
            )
        return
    if isinstance(node, ast.Name):
        if id(node) in callees:
            return  # checked as the callee of its Call
        if node.id not in allowed:
            raise ExpressionError(
                f"{where}: unknown name {node.id!r} in {text!r}; "
                f"known: {sorted(allowed)}"
            )
        used.add(node.id)
        return
    raise ExpressionError(f"{where}: {type(node).__name__} is not allowed in {text!r}")
