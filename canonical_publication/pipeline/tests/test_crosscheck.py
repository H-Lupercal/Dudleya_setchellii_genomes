import math

import pytest
from organelle_pipeline.crosscheck import ratio_of_jointly_defined_components


def test_hudson_crosscheck_excludes_denominator_when_numerator_is_undefined() -> None:
    result = ratio_of_jointly_defined_components(
        (math.nan, 2.0, 0.0),
        (9.2, 4.0, 0.0),
    )

    assert result.numerator == 2.0
    assert result.denominator == 4.0
    assert result.jointly_defined_sites == 2
    assert result.ratio == 0.5


def test_hudson_crosscheck_rejects_misaligned_component_arrays() -> None:
    with pytest.raises(ValueError, match="equal lengths"):
        ratio_of_jointly_defined_components((1.0,), (1.0, 2.0))
