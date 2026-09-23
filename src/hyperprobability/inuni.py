"""In-universe probability: ordinary probability inside a contained universe.

These are the probabilities a run can see from inside the universe: the
chance of an event at a given stage, the chance of meeting it within the
run's own finite time, and the long-run frequency with which it holds. They
are real numbers in [0, 1], exact rationals here. Hyperprobability measures
something else: how many evaluations make the event certain.
"""

from __future__ import annotations

from fractions import Fraction

from ordinatics.ordinals import OMEGA

from hyperprobability.transfinite import TransfiniteLaw
from hyperprobability.universe import Universe, _stage


def probability(universe: Universe, event: object, initial: object, stage: object = 0) -> Fraction:
    """``P_μ(X_stage ∈ E)``, the chance that the run is in ``event`` at ``stage``.

    At a finite stage this is ordinary probability. At a transfinite stage it
    is transfinite probability, from :meth:`Universe.law`.
    """
    marked = universe._event(event)
    law = universe._push(universe._vector(initial), _stage(stage))
    return sum((p for q, p in law.items() if q in marked), Fraction(0))


def reach_probability(universe: Universe, event: object, initial: object) -> Fraction:
    """``P_μ(τ_E < ω)``: the chance that the run meets ``event`` at some finite stage."""
    return TransfiniteLaw(universe, event, initial).before(OMEGA)


def frequency(universe: Universe, event: object, initial: object) -> Fraction:
    """The expected long-run fraction of finite stages at which the run is in ``event``.

    Almost surely the run settles in an attractor ``A`` and then spends the
    fraction ``π_A(E)`` of its time in ``E``, by the ergodic theorem. This is
    the expectation ``Σ_A P(settle in A)·π_A(E)``, the mass of ``event`` under
    :meth:`Universe.in_uni`.
    """
    marked = universe._event(event)
    law = universe._in_uni(universe._vector(initial))
    return sum((p for q, p in law.items() if q in marked), Fraction(0))
