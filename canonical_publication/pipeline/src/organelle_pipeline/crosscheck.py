"""Small, implementation-independent helpers for trusted statistic checks."""

from __future__ import annotations

import math
from collections.abc import Iterable
from dataclasses import dataclass


@dataclass(frozen=True)
class JointComponentRatio:
    numerator: float
    denominator: float
    jointly_defined_sites: int
    ratio: float


def ratio_of_jointly_defined_components(
    numerators: Iterable[float],
    denominators: Iterable[float],
) -> JointComponentRatio:
    """Sum only sites where both ratio components are finite.

    Hudson numerators require at least two allele calls within each
    population.  A between-population denominator can still be finite with
    only one call, but including that denominator when its numerator is
    undefined biases the ratio downward.
    """

    numerator_values = list(numerators)
    denominator_values = list(denominators)
    if len(numerator_values) != len(denominator_values):
        raise ValueError("Numerator and denominator arrays must have equal lengths")
    pairs = [
        (numerator, denominator)
        for numerator, denominator in zip(numerator_values, denominator_values, strict=True)
        if math.isfinite(numerator) and math.isfinite(denominator)
    ]
    numerator = sum(value[0] for value in pairs)
    denominator = sum(value[1] for value in pairs)
    return JointComponentRatio(
        numerator=numerator,
        denominator=denominator,
        jointly_defined_sites=len(pairs),
        ratio=numerator / denominator if denominator else math.nan,
    )
