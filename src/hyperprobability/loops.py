"""Strange loops, the Kac bridge between the two probabilities, and self-simulation.

A *strange loop* is an attractor whose limit does not collapse to its own
in-universe law: the level above sends a run that settles there somewhere
else, perhaps back to where it began. Strange loops are the only way a
transfinite law can jump at a limit stage (SPEC Theorem 5.2).

The *Kac bridge* relates in-universe probability to hyperprobability inside
a level-0 attractor ``A`` that meets an event ``E``:

    π_A(E) ≥ 1/ι(ℍ⁺_A(E)),

where ``π_A(E)`` is the long-run fraction of time the run spends in ``E`` and
``ℍ⁺_A(E)`` is the worst guarantee of a return to ``E``. So the ordinal
probability of a start spread over ``A ∩ E`` is a lower bound on in-universe
probability, sharp exactly when every return takes the same time (SPEC
Theorem 6.2). A single start can exceed it (SPEC Proposition 7.5).

:func:`refine` is the guarantee of a simulation carried out by a finer one
(SPEC Proposition 8.1), and :func:`self_simulation` the guarantee of a
simulation that has to run itself (SPEC Proposition 8.3).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from fractions import Fraction
from typing import Optional

from ordinatics.calculus import least_fixed_point
from ordinatics.ordinals import ZERO, Ordinal

from hyperprobability.hyper import INFINITY, Guarantee, _return_guarantee
from hyperprobability.ordinal import OmegaFraction
from hyperprobability.universe import Attractor, Universe, _count, _stage


def strange_loops(universe: Universe, level: Optional[int] = None) -> tuple[Attractor, ...]:
    """The attractors whose limit is not their own in-universe law.

    With ``level``, only that level's attractors are searched. Only the
    levels below ``universe.depth`` have limit rules, so only they can loop.
    """
    levels = range(universe.depth) if level is None else (_count(level, "level"),)
    return tuple(a for k in levels for a in universe.attractors(k) if not a.collapses)


@dataclass(frozen=True)
class KacBound:
    """The Kac bridge for one level-0 attractor that meets an event.

    ``event`` holds the states of the event inside ``attractor``.
    ``frequency`` is the in-universe probability ``π_A(E)``, the long-run
    fraction of time in the event. ``guarantee`` is ``ℍ⁺_A(E)``, the largest
    return hyperprobability from those states, and ``ordinal_probability`` is
    ``1/ι(ℍ⁺_A(E))``.
    """

    attractor: Attractor
    event: frozenset
    frequency: Fraction
    guarantee: Ordinal
    ordinal_probability: OmegaFraction

    @property
    def holds(self) -> bool:
        """Whether ``ordinal_probability ≤ frequency``, which SPEC Theorem 6.2 proves always holds."""
        return self.ordinal_probability <= self.frequency

    @property
    def tight(self) -> bool:
        """Whether the bound is an equality: every return takes exactly ``guarantee`` steps."""
        return self.ordinal_probability == self.frequency


def kac(universe: Universe, event: object) -> tuple[KacBound, ...]:
    """The Kac bridge for each level-0 attractor that meets ``event``.

    Inside such an attractor ``A`` the run returns to ``E`` again and again.
    By Kac's lemma its mean return time, from the stationary law on
    ``A ∩ E``, is ``1/π_A(E)``. No return takes longer than ``ℍ⁺_A(E)``, so
    ``π_A(E) ≥ 1/ι(ℍ⁺_A(E))``. When returns are certain but unbounded,
    ``ℍ⁺_A(E) = ω`` and the bound is the infinitesimal ``1/ω``.
    """
    marked = universe._event(event)
    data = universe._level(0)
    bounds = []
    for attractor, members, law in zip(universe.attractors(0), data.classes, data.stationary):
        inside = {q: law[q] for q in members if q in marked}
        if not inside:
            continue
        guarantee = _return_guarantee(universe, marked, inside)
        if guarantee is INFINITY:  # pragma: no cover - a recurrent class always returns
            raise AssertionError(f"no guaranteed return inside {attractor!r}")
        bounds.append(
            KacBound(
                attractor=attractor,
                event=universe._members(sorted(inside)),
                frequency=sum(inside.values(), Fraction(0)),
                guarantee=guarantee,
                ordinal_probability=OmegaFraction.from_ordinal(guarantee).reciprocal(),
            )
        )
    return tuple(bounds)


def refine(guarantee: Guarantee, cost: object) -> Guarantee:
    """The guarantee of a run carried out at ``cost`` evaluations per step: ``cost·guarantee``.

    If each stage of a run takes ``cost`` stages of a finer simulation, run
    stage ``α`` is simulation stage ``cost·α``, the ordinal product: it adds
    ``cost`` at each successor and is continuous at limits. So a finite cost
    leaves an infinite guarantee alone, since ``2·ω = ω``, while a cost of
    ``ω`` doubles ``2`` into ``ω·2``. ``guarantee`` is an Ordinal, a natural
    number or :data:`INFINITY`. For a finite cost this is the guarantee of the
    subdivided universe (SPEC Proposition 8.1); an infinite cost is an
    interpretation (SPEC Remark 8.2).
    """
    factor = _stage(cost)
    if factor == ZERO:
        raise ValueError("a simulation step costs at least one evaluation")
    if guarantee is INFINITY:
        return INFINITY
    try:
        alpha = _stage(guarantee)
    except TypeError:
        raise TypeError("a guarantee is an Ordinal, a natural number or INFINITY") from None
    return factor * alpha


def self_simulation(cost: Callable[[Ordinal], object], start: object = 0) -> Ordinal:
    """The least guarantee ``α ≥ start`` that pays for itself: ``cost(α) == α``.

    A run that must first simulate itself needs a guarantee ``α`` covering
    that simulation, so ``α = cost(α)``. For ``cost(α) = h + α``, meaning
    ``h`` evaluations and then the whole run again, the answer is ``h·ω``:
    ``ω`` when ``h = 1``, and ``ω**2`` when ``h = ω``. ``cost`` must be weakly
    increasing. The search is :func:`ordinatics.calculus.least_fixed_point`,
    whose result is a checked fixed point but whose passage through limits
    is heuristic. An interpretation (SPEC Proposition 8.3).
    """
    return least_fixed_point(lambda alpha: _stage(cost(alpha)), _stage(start))
