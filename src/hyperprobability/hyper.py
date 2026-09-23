"""Hyperprobability: the fewest simulation evaluations that guarantee an event.

A *simulation evaluation* is one stage of a run of a universe: a step of its
kernel, or a limit stage. It is not hypermath's simulation relation ``==``.
For an event ``E`` and an initial law ``μ``, the hyperprobability is

    ℍ_μ(E) = the least ordinal α with P_μ(τ_E ≤ α) = 1,

where ``τ_E`` is the first stage at which the run is in ``E``. It is ``∞``
when no stage below ``ω**ω`` guarantees ``E``. Guarantees come in transfinite
kinds, which is why the measure is an ordinal:

* a clock guarantees its next tick within 2 evaluations;
* a fair coin guarantees heads only at ``ω``, which is certain by then and
  never by any finite stage;
* a strange loop can guarantee an event at a limit stage that no finite run
  ever reaches.

See SPEC.md §4.
"""

from __future__ import annotations

from collections.abc import Hashable
from typing import Union

from ordinatics.ordinals import ONE, ZERO, Ordinal

from hyperprobability._exact import Row, reachable
from hyperprobability.universe import Universe, _Level


class Infinity:
    """The hyperprobability of an event that no stage guarantees.

    It is above every ordinal and absorbs ordinal addition. There is one
    instance, :data:`INFINITY`.
    """

    __slots__ = ()
    _instance: Union[Infinity, None] = None

    def __new__(cls) -> Infinity:
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    @staticmethod
    def _comparable(other: object) -> bool:
        if isinstance(other, bool):
            return False
        return isinstance(other, (Ordinal, Infinity, int))

    def __eq__(self, other: object) -> bool:
        return other is self

    def __hash__(self) -> int:
        return hash("hyperprobability.INFINITY")

    def __lt__(self, other: object) -> bool:
        if not self._comparable(other):
            return NotImplemented
        return False

    def __le__(self, other: object) -> bool:
        if not self._comparable(other):
            return NotImplemented
        return other is self

    def __gt__(self, other: object) -> bool:
        if not self._comparable(other):
            return NotImplemented
        return other is not self

    def __ge__(self, other: object) -> bool:
        if not self._comparable(other):
            return NotImplemented
        return True

    def __add__(self, other: object) -> Infinity:
        if not self._comparable(other):
            return NotImplemented
        return self

    __radd__ = __add__

    def __reduce__(self) -> str:
        return "INFINITY"

    def __repr__(self) -> str:
        return "INFINITY"

    def __str__(self) -> str:
        return "∞"


INFINITY = Infinity()

Guarantee = Union[Ordinal, Infinity]


def hyperprobabilities(universe: Universe, event: object) -> dict[Hashable, Guarantee]:
    """``ℍ_q(E)`` for every state ``q`` of ``universe``."""
    return dict(zip(universe.states, _guarantees(universe, universe._event(event))))


def hyperprobability(universe: Universe, event: object, initial: object) -> Guarantee:
    """``ℍ_μ(E)``: the least stage by which ``event`` has surely happened.

    ``initial`` is a state or a law on states. ``ℍ_μ(E)`` is the largest
    ``ℍ_q(E)`` over the states ``q`` that ``μ`` can start in (SPEC
    Proposition 4.3).
    """
    values = _guarantees(universe, universe._event(event))
    return max(values[q] for q in universe._vector(initial))


def return_hyperprobability(universe: Universe, event: object, initial: object) -> Guarantee:
    """``ℍ⁺_μ(E)``: the least stage ``α ≥ 1`` by which the run has surely been in ``event``.

    Unlike ``ℍ``, the start itself does not count, so for a start inside the
    event this is the guarantee of a return. It equals ``1 + max ℍ_r(E)`` over
    the states ``r`` one step from the start, with ordinal addition, so
    ``1 + ω == ω`` (SPEC Proposition 4.8).
    """
    return _return_guarantee(universe, universe._event(event), universe._vector(initial))


def _return_guarantee(universe: Universe, marked: frozenset[int], vector: Row) -> Guarantee:
    values = _guarantees(universe, marked)
    rows = universe._level(0).rows
    worst = max(values[r] for q in vector for r in rows[q])
    return worst if worst is INFINITY else ONE + worst


def _guarantees(universe: Universe, marked: frozenset[int]) -> list[Guarantee]:
    """``ℍ_q(E)`` for every state index ``q``, by the level recursion of SPEC Theorem 4.4."""
    key = ("guarantees", marked)
    values = universe._memo.get(key)
    if values is None:
        stopped = universe._stop_at(marked)
        count = len(universe.states)
        values = [ZERO if q in marked else INFINITY for q in range(count)]
        for k in range(universe.depth + 1):
            level = stopped._level(k)
            values = [_resolve(q, k, level, marked, values) for q in range(count)]
        universe._memo[key] = values
    return values


def _resolve(
    q: int, k: int, level: _Level, marked: frozenset[int], previous: list[Guarantee]
) -> Guarantee:
    """``R_k(q)`` from ``R_(k-1)``: the least stage below ``ω**(k+1)`` that guarantees the event.

    Stage ``ω**k·c + β`` guarantees it from ``q`` exactly when every state the
    level-``k`` kernel can reach from ``q`` in ``c`` steps is in the event or
    is guaranteed from within ``β < ω**k``. Once some ``c`` works every larger
    one does, and then one below the number of states does, so the search is
    short. If no stage below ``ω**(k+1)`` works, stage ``ω**(k+1)`` works
    exactly when the level-``k`` run surely meets the event.
    """
    bound = Ordinal.omega_power(k)
    support = frozenset((q,))
    for count in range(len(level.succ)):
        pending = [previous[r] for r in support if r not in marked]
        if all(value < bound for value in pending):
            return bound * count + max(pending, default=ZERO)
        support = frozenset(t for r in support for t in level.succ[r])
    reach = reachable((q,), level.succ)
    for members in level.classes:
        if members[0] in reach and members[0] not in marked:
            return INFINITY
    return Ordinal.omega_power(k + 1)
