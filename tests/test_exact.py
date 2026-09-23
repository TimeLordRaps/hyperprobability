import random
from fractions import Fraction

import pytest

from hyperprobability._exact import (
    absorption,
    closed_classes,
    compose,
    period,
    push,
    push_power,
    reachable,
    solve,
    stationary,
    strongly_connected,
)


def random_rows(seed: int, n: int) -> list[dict[int, Fraction]]:
    rng = random.Random(seed)
    rows = []
    for _ in range(n):
        chosen = rng.sample(range(n), rng.randint(1, min(n, 3)))
        weights = [rng.randint(1, 5) for _ in chosen]
        rows.append({t: Fraction(w, sum(weights)) for t, w in zip(chosen, weights)})
    return rows


def successors(rows):
    return [sorted(row) for row in rows]


def reach_matrix(n, succ):
    return [reachable((q,), succ) for q in range(n)]


def test_solve_returns_the_exact_solution():
    rng = random.Random(1)
    solved = 0
    while solved < 20:
        n = rng.randint(1, 5)
        matrix = [[Fraction(rng.randint(-4, 4)) for _ in range(n)] for _ in range(n)]
        x = [[Fraction(rng.randint(-5, 5), rng.randint(1, 5))] for _ in range(n)]
        rhs = [[sum(matrix[i][k] * x[k][0] for k in range(n))] for i in range(n)]
        try:
            result = solve(matrix, rhs)
        except ArithmeticError:
            continue
        assert result == x
        solved += 1


def test_solve_rejects_a_singular_system():
    with pytest.raises(ArithmeticError):
        solve(
            [[Fraction(1), Fraction(2)], [Fraction(2), Fraction(4)]], [[Fraction(1)], [Fraction(2)]]
        )


def test_components_and_closed_classes_match_reachability():
    for seed in range(40):
        n = random.Random(seed).randint(1, 7)
        succ = successors(random_rows(seed, n))
        reach = reach_matrix(n, succ)
        components = strongly_connected(n, succ)
        assert sorted(q for c in components for q in c) == list(range(n))
        for component in components:
            for q in component:
                for r in range(n):
                    assert (r in component) == (r in reach[q] and q in reach[r])
        closed = closed_classes(n, succ)
        expected = [c for c in components if all(reach[q] <= set(c) for q in c)]
        assert sorted(map(sorted, closed)) == sorted(map(sorted, expected))
        assert [min(c) for c in closed] == sorted(min(c) for c in closed)


def test_period():
    cycle = [[1], [2], [0]]
    assert period([0, 1, 2], cycle) == 3
    assert period([0, 1, 2], [[1], [2], [0, 2]]) == 1
    assert period([0, 1, 2, 3], [[1], [0, 2], [3], [2, 0]]) == 2
    assert period([0], [[0]]) == 1


def test_stationary_laws_are_invariant():
    for seed in range(40):
        n = random.Random(seed).randint(1, 6)
        rows = random_rows(seed, n)
        for members in closed_classes(n, successors(rows)):
            law = stationary(rows, members)
            assert sum(law.values()) == 1
            assert set(law) == set(members)
            assert push(law, rows) == law


def test_absorption_is_harmonic_and_complete():
    for seed in range(40):
        n = random.Random(seed).randint(1, 6)
        rows = random_rows(seed, n)
        classes = closed_classes(n, successors(rows))
        h = absorption(rows, classes)
        for q in range(n):
            assert sum(h[q]) == 1
            for a, members in enumerate(classes):
                assert h[q][a] == sum(p * h[t][a] for t, p in rows[q].items())
                if q in members:
                    assert h[q][a] == 1


def test_push_power_matches_repeated_steps():
    rows = random_rows(7, 4)
    vector = {0: Fraction(1, 3), 2: Fraction(2, 3)}
    expected = dict(vector)
    for times in range(70):
        assert push_power(vector, rows, times) == expected
        expected = push(expected, rows)


def test_compose_is_the_kernel_product():
    first, second = random_rows(3, 4), random_rows(4, 4)
    product = compose(first, second)
    for q in range(4):
        assert product[q] == push(push({q: Fraction(1)}, first), second)


@pytest.mark.parametrize("n", [0, 1])
def test_trivial_sizes(n):
    assert strongly_connected(n, [[i] for i in range(n)]) == [[i] for i in range(n)]
