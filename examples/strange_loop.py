"""Strange loops: guarantees at limit stages that no finite run reaches.

Each universe starts at ``"a"`` and asks for the event ``"out"``, which no
finite run reaches from there. In-universe probability gives it 0 every time.
Hyperprobability tells the four universes apart:

* ``dead end``: ``"a"`` is a fixed point whose limit collapses, so ``out``
  never happens: ``ℍ = ∞``;
* ``loop``: a two-cycle whose limit rule leaves through ``out``: ``ℍ = ω``;
* ``chain``: the same loop, whose limit starts a three-step walk: ``ℍ = ω + 3``;
* ``swap``: two fixed points that swap at every limit, until the limit of
  the swapping leaves: ``ℍ = ω^2``.

Run it with ``python examples/strange_loop.py``.
"""

from hyperprobability import (
    INFINITY,
    TransfiniteLaw,
    Universe,
    hyperprobability,
    ordinal_probability,
    reach_probability,
    strange_loops,
)

pair = frozenset({"a", "b"})
universes = {
    "dead end": Universe.deterministic({"a": "a", "out": "out"}),
    "loop": Universe.deterministic({"a": "b", "b": "a", "out": "out"}, [{pair: "out"}]),
    "chain": Universe.deterministic(
        {"a": "b", "b": "a", "x1": "x2", "x2": "x3", "x3": "out", "out": "out"}, [{pair: "x1"}]
    ),
    "swap": Universe.deterministic(
        {"a": "a", "b": "b", "out": "out"},
        [{frozenset({"a"}): "b", frozenset({"b"}): "a"}, {pair: "out"}],
    ),
}


def show(label: str, value: object) -> None:
    print(f"  {label:<30}{value}")


for name, universe in universes.items():
    law = TransfiniteLaw(universe, "out", "a")
    guarantee = hyperprobability(universe, "out", "a")
    loops = [f"{sorted(a.states)} at level {a.level}" for a in strange_loops(universe)]
    print(f"{name}:")
    show("strange loops", ", ".join(loops) or "none")
    show("reached in finite time", reach_probability(universe, "out", "a"))
    show("hyperprobability", guarantee)
    show("ordinal probability", ordinal_probability(universe, "out", "a"))
    show("never reached", law.never())
    if guarantee is not INFINITY:
        show(f"first reached at stage {guarantee}", law.exactly(guarantee))
