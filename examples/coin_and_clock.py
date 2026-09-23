"""A clock and a coin: the same in-universe probability, different guarantees.

Both universes spend half their time in the event, so in-universe probability
cannot tell them apart. Hyperprobability can: the clock surely returns to its
event within 2 evaluations, while the coin surely returns only by stage ``ω``.
Ordinal probability ``1/ι(ℍ⁺)`` is then ``1/2`` for the clock and the
infinitesimal ``1/ω`` for the coin, and the Kac bridge bounds each frequency
from below.

Run it with ``python examples/coin_and_clock.py``.
"""

from hyperprobability import (
    OMEGA,
    TransfiniteLaw,
    Universe,
    frequency,
    hyperprobability,
    kac,
    ordinal_probability,
    return_hyperprobability,
)

clock = Universe.deterministic({"tick": "tock", "tock": "tick"})
coin = Universe({"H": {"H": "1/2", "T": "1/2"}, "T": {"H": "1/2", "T": "1/2"}})

for name, universe, event, other in [("clock", clock, "tick", "tock"), ("coin", coin, "H", "T")]:
    (bound,) = kac(universe, event)
    print(f"{name}, event {event!r}:")
    print(f"  in-universe frequency       {frequency(universe, event, other)}")
    print(f"  hyperprobability from {other!r}  {hyperprobability(universe, event, other)}")
    print(f"  return hyperprobability     {return_hyperprobability(universe, event, event)}")
    print(f"  ordinal probability         {ordinal_probability(universe, event, event)}")
    print(
        f"  Kac bridge                  {bound.ordinal_probability} <= {bound.frequency}, "
        f"tight: {bound.tight}"
    )

law = TransfiniteLaw(coin, "H", "T")
first = ", ".join(str(law.exactly(n)) for n in range(1, 6))
print(f"coin: first heads at stage 1, 2, 3, 4, 5 with probability {first}, ...")
print(f"coin: heads by stage ω with probability {law.at_most(OMEGA)}")
