"""Exact rational linear algebra and graph routines for finite kernels.

Every number here is a :class:`fractions.Fraction`; nothing is rounded. A
kernel is a list of rows, one per state index, and a row maps a successor's
index to its probability. Rows keep only nonzero entries.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from fractions import Fraction
from math import gcd

Row = dict[int, Fraction]

_ZERO = Fraction(0)
_ONE = Fraction(1)


def solve(matrix: list[list[Fraction]], rhs: list[list[Fraction]]) -> list[list[Fraction]]:
    """Return X with ``matrix @ X == rhs`` for a nonsingular square ``matrix``.

    Gauss-Jordan elimination over the rationals. ``rhs`` has one row per
    equation and any number of columns.
    """
    n = len(matrix)
    table = [list(matrix[i]) + list(rhs[i]) for i in range(n)]
    for col in range(n):
        pivot = next((r for r in range(col, n) if table[r][col]), None)
        if pivot is None:
            raise ArithmeticError("singular linear system")
        table[col], table[pivot] = table[pivot], table[col]
        head = table[col][col]
        if head != 1:
            table[col] = [x / head for x in table[col]]
        for r in range(n):
            factor = table[r][col]
            if r != col and factor:
                table[r] = [x - factor * y for x, y in zip(table[r], table[col])]
    return [row[n:] for row in table]


def strongly_connected(n: int, succ: Sequence[Sequence[int]]) -> list[list[int]]:
    """Strongly connected components of a digraph on ``range(n)`` (iterative Tarjan)."""
    index: list[int | None] = [None] * n
    low = [0] * n
    on_stack = [False] * n
    stack: list[int] = []
    components: list[list[int]] = []
    counter = 0
    for root in range(n):
        if index[root] is not None:
            continue
        index[root] = low[root] = counter
        counter += 1
        stack.append(root)
        on_stack[root] = True
        work = [(root, 0)]
        while work:
            v, i = work[-1]
            if i < len(succ[v]):
                work[-1] = (v, i + 1)
                w = succ[v][i]
                if index[w] is None:
                    index[w] = low[w] = counter
                    counter += 1
                    stack.append(w)
                    on_stack[w] = True
                    work.append((w, 0))
                elif on_stack[w]:
                    low[v] = min(low[v], index[w])  # type: ignore[type-var]
                continue
            work.pop()
            if work:
                parent = work[-1][0]
                low[parent] = min(low[parent], low[v])
            if low[v] == index[v]:
                component = []
                while True:
                    w = stack.pop()
                    on_stack[w] = False
                    component.append(w)
                    if w == v:
                        break
                components.append(sorted(component))
    return components


def closed_classes(n: int, succ: Sequence[Sequence[int]]) -> list[list[int]]:
    """The closed strongly connected components, ordered by their least index."""
    result = []
    for component in strongly_connected(n, succ):
        members = set(component)
        if all(t in members for s in component for t in succ[s]):
            result.append(component)
    result.sort(key=lambda c: c[0])
    return result


def period(members: Sequence[int], succ: Sequence[Sequence[int]]) -> int:
    """The period of a strongly connected class: the gcd of its cycle lengths."""
    inside = set(members)
    level = {members[0]: 0}
    queue = [members[0]]
    result = 0
    for v in queue:
        for w in succ[v]:
            if w not in inside:
                continue
            if w in level:
                result = gcd(result, level[v] + 1 - level[w])
            else:
                level[w] = level[v] + 1
                queue.append(w)
    return result


def reachable(sources: Iterable[int], succ: Sequence[Sequence[int]]) -> set[int]:
    """Every index reachable from ``sources`` along edges, sources included."""
    seen = set(sources)
    queue = list(seen)
    for v in queue:
        for w in succ[v]:
            if w not in seen:
                seen.add(w)
                queue.append(w)
    return seen


def stationary(rows: Sequence[Row], members: Sequence[int]) -> Row:
    """The unique stationary law of an irreducible closed class."""
    k = len(members)
    if k == 1:
        return {members[0]: _ONE}
    position = {s: i for i, s in enumerate(members)}
    # pi (P - I) = 0 is (P^T - I) pi^T = 0; the last equation becomes sum(pi) = 1.
    matrix = [[_ZERO] * k for _ in range(k)]
    for s in members:
        column = position[s]
        for t, p in rows[s].items():
            matrix[position[t]][column] += p
    for i in range(k):
        matrix[i][i] -= 1
    matrix[k - 1] = [_ONE] * k
    rhs = [[_ZERO] for _ in range(k)]
    rhs[k - 1][0] = _ONE
    solution = solve(matrix, rhs)
    return {members[i]: solution[i][0] for i in range(k)}


def absorption(rows: Sequence[Row], classes: Sequence[Sequence[int]]) -> list[list[Fraction]]:
    """``h[q][a]``: the probability that the chain from ``q`` ends in ``classes[a]``.

    In a finite chain every run is eventually absorbed in some closed class, so
    each row sums to one.
    """
    n = len(rows)
    width = len(classes)
    h = [[_ZERO] * width for _ in range(n)]
    home: dict[int, int] = {}
    for a, members in enumerate(classes):
        for s in members:
            home[s] = a
            h[s][a] = _ONE
    transient = [q for q in range(n) if q not in home]
    if not transient:
        return h
    position = {q: i for i, q in enumerate(transient)}
    t = len(transient)
    matrix = [[_ZERO] * t for _ in range(t)]
    rhs = [[_ZERO] * width for _ in range(t)]
    for q in transient:
        i = position[q]
        matrix[i][i] += 1
        for r, p in rows[q].items():
            if r in position:
                matrix[i][position[r]] -= p
            else:
                rhs[i][home[r]] += p
    solution = solve(matrix, rhs)
    for q in transient:
        h[q] = solution[position[q]]
    return h


def push(vector: Row, rows: Sequence[Row]) -> Row:
    """One step of a law: ``vector @ rows``."""
    result: Row = {}
    for j, p in vector.items():
        for k, q in rows[j].items():
            result[k] = result.get(k, _ZERO) + p * q
    return {k: v for k, v in result.items() if v}


def compose(first: Sequence[Row], second: Sequence[Row]) -> list[Row]:
    """The kernel product ``first @ second``: first ``first``, then ``second``."""
    return [push(row, second) for row in first]


def push_power(vector: Row, rows: Sequence[Row], times: int) -> Row:
    """``vector @ rows ** times``, by repeated squaring once ``times`` is large."""
    if times <= 32:
        for _ in range(times):
            vector = push(vector, rows)
        return vector
    square = list(rows)
    while times:
        if times & 1:
            vector = push(vector, square)
        times >>= 1
        if times:
            square = compose(square, square)
    return vector
