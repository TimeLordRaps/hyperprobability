"""Shared universes and brute-force oracles for the tests.

The oracles search stages directly with exact laws and share no code with the
level recursion of SPEC Theorem 4.4, which is what they check.
"""

from __future__ import annotations

import itertools
import random
from fractions import Fraction

from hyperprobability import INFINITY, ONE, Ordinal, Universe

HALF = Fraction(1, 2)


def clock() -> Universe:
    """Two states that alternate forever."""
    return Universe.deterministic({"tick": "tock", "tock": "tick"})


def coin() -> Universe:
    """A fair coin tossed at every step."""
    toss = {"H": HALF, "T": HALF}
    return Universe({"H": toss, "T": toss})


def loop() -> Universe:
    """A two-cycle whose limit is a way out that no finite run reaches."""
    return Universe.deterministic(
        {"a": "b", "b": "a", "out": "out"}, [{frozenset({"a", "b"}): "out"}]
    )


def chain() -> Universe:
    """The same loop, whose limit starts a three-step path out."""
    return Universe.deterministic(
        {"a": "b", "b": "a", "x1": "x2", "x2": "x3", "x3": "out", "out": "out"},
        [{frozenset({"a", "b"}): "x1"}],
    )


def swap() -> Universe:
    """Two fixed points that swap at every limit, until the limit of the swapping exits."""
    return Universe.deterministic(
        {"a": "a", "b": "b", "out": "out"},
        [
            {frozenset({"a"}): "b", frozenset({"b"}): "a"},
            {frozenset({"a", "b"}): "out"},
        ],
    )


def fork() -> Universe:
    """The same loop, whose limit forks onto a one-step and a two-step path out."""
    return Universe.deterministic(
        {"a": "b", "b": "a", "x": "out", "y1": "y2", "y2": "out", "out": "out"},
        [{frozenset({"a", "b"}): {"x": HALF, "y1": HALF}}],
    )


def dead_end() -> Universe:
    """A state that never leaves, beside the event it never meets."""
    return Universe.deterministic({"a": "a", "out": "out"})


def random_law(rng: random.Random, states: list) -> dict:
    chosen = rng.sample(states, rng.randint(1, len(states)))
    weights = [rng.randint(1, 4) for _ in chosen]
    total = sum(weights)
    return {state: Fraction(weight, total) for state, weight in zip(chosen, weights)}


def random_rule(seed: str, states: list):
    """A limit rule that sends each attractor somewhere fixed by ``seed``."""

    def rule(members: frozenset):
        rng = random.Random(f"{seed}/{sorted(members)}")
        roll = rng.random()
        if roll < 0.2:
            return None
        if roll < 0.7:
            return rng.choice(states)
        return random_law(rng, states)

    return rule


def random_universe(seed: int, size: int = 0, depth: int = -1) -> Universe:
    """A small random universe with absorbing, deterministic and random rows.

    Most rows and limit targets are single states, so that events are often
    met late: after a strange loop, on a walk that follows one, or at the
    limit of a loop of loops. Half the universes are almost deterministic and
    a little larger, which is where stages like ``ω·2 + 3`` come from.
    """
    rng = random.Random(seed)
    walks = rng.random() < 0.5
    states = list(range(size or rng.randint(1, 6 if walks else 4)))
    kernel = {}
    for state in states:
        roll = rng.random()
        if roll < 0.12:
            kernel[state] = {state: 1}
        elif roll < (0.92 if walks else 0.6):
            kernel[state] = {rng.choice(states): 1}
        else:
            kernel[state] = random_law(rng, states)
    levels = rng.randint(0, 3) if depth < 0 else depth
    return Universe(kernel, [random_rule(f"{seed}/{k}", states) for k in range(levels)])


def random_events(universe: Universe, seed: int, count: int = 3) -> list[frozenset]:
    """The empty event, the whole universe, and a few small random events."""
    rng = random.Random(seed)
    states = list(universe.states)
    events = [frozenset(), frozenset(states)]
    for _ in range(count):
        events.append(frozenset(rng.sample(states, rng.randint(1, max(1, len(states) // 2)))))
    return events


def ordinals(exponents: int, bound: int):
    """Every ordinal below ``ω**exponents`` with coefficients at most ``bound``, in order."""
    for digits in itertools.product(range(bound + 1), repeat=exponents):
        yield Ordinal(tuple(reversed(digits)))


def mass(law: dict, event: frozenset) -> Fraction:
    return sum((p for state, p in law.items() if state in event), Fraction(0))


def first_certain(universe: Universe, event: frozenset, initial: object):
    """The least stage ``α`` with ``P(τ_E ≤ α) = 1``, by search over stages in order.

    ``P(τ_E ≤ α)`` is the mass of the event under the law at ``α`` of the
    stopped universe, and it never decreases in ``α``. So a block of stages
    whose supremum is not certain holds no certain stage and is skipped.
    Below ``ω**(depth+1)`` the least certain stage has coefficients below
    ``2**n``: certainty depends only on which states the law can be in, and
    those sets repeat within ``2**n`` steps at any level.
    """
    stopped = universe.stop_at(event)
    top = universe.depth + 1
    bound = 2 ** len(universe) + 1

    def search(law: dict, exponent: int, prefix: Ordinal):
        if exponent < 0:
            return prefix if mass(law, event) == 1 else None
        step = Ordinal.omega_power(exponent)
        for digit in range(bound + 1):
            after = stopped.law(law, step)
            if mass(after, event) == 1:
                found = search(law, exponent - 1, prefix + step * digit)
                if found is not None:
                    return found
            law = after
        return None

    limit = Ordinal.omega_power(top)
    if mass(stopped.law(initial, limit), event) != 1:
        return INFINITY
    found = search(stopped.law(initial), top - 1, Ordinal.from_int(0))
    return limit if found is None else found


def first_certain_return(universe: Universe, event: frozenset, state: object):
    """The least stage ``α ≥ 1`` by which the run from ``state`` has surely been in ``event``."""
    later = first_certain(universe, event, universe.step(state))
    return later if later is INFINITY else ONE + later
