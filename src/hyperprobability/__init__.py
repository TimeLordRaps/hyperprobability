"""Hyperprobability: the fewest simulation evaluations that guarantee an event.

Probability is in-universe: inside a contained universe it is the ordinary
chance of an event. Hyperprobability ``ℍ_μ(E)`` measures something else: the
least ordinal stage by which a run of the universe has surely produced the
event, counted in transfinite evaluations. Transfinite probability, the law
of the first stage at which the event happens, collapses out through the
universe's attractors, and jumps only through strange loops. Ordinal
probability ``1/ι(ℍ⁺)`` in the ordered field ℚ(ω) is built on top of
hyperprobability.

Ordinals come from :mod:`ordinatics`. See SPEC.md for definitions, proofs
and the status of every claim.
"""

from ordinatics.ordinals import OMEGA, ONE, ZERO, Ordinal

from hyperprobability.hyper import (
    INFINITY,
    Guarantee,
    Infinity,
    hyperprobabilities,
    hyperprobability,
    return_hyperprobability,
)
from hyperprobability.inuni import frequency, probability, reach_probability
from hyperprobability.loops import KacBound, kac, refine, self_simulation, strange_loops
from hyperprobability.ordinal import OmegaFraction, ordinal_probability
from hyperprobability.transfinite import TransfiniteLaw
from hyperprobability.universe import Attractor, Universe, exact

__version__ = "0.1.0"

__all__ = [
    "INFINITY",
    "OMEGA",
    "ONE",
    "ZERO",
    "Attractor",
    "Guarantee",
    "Infinity",
    "KacBound",
    "OmegaFraction",
    "Ordinal",
    "TransfiniteLaw",
    "Universe",
    "exact",
    "frequency",
    "hyperprobabilities",
    "hyperprobability",
    "kac",
    "ordinal_probability",
    "probability",
    "reach_probability",
    "refine",
    "return_hyperprobability",
    "self_simulation",
    "strange_loops",
]
