import copy
import tomllib
from pathlib import Path

import pytest
from organelle_pipeline.configuration import (
    PublicationConfigurationError,
    validate_publication_config,
)


def canonical_config() -> dict[str, object]:
    root = Path(__file__).resolve().parents[3]
    return tomllib.loads((root / "canonical_publication/config/publication_config.toml").read_text())


def test_current_publication_configuration_is_accepted() -> None:
    validate_publication_config(canonical_config())


@pytest.mark.parametrize(
    ("section", "key", "unsafe_value"),
    [
        ("mapping", "minimum_mapping_quality", 10),
        ("variants", "minimum_genotype_quality", 0),
        ("variants", "primary_minimum_minor_allele_count", 2),
        ("phylogeny", "ultrafast_bootstrap_replicates", 100),
        ("admixture", "maximum_k", 8),
    ],
)
def test_acceptance_critical_policy_cannot_be_silently_weakened(
    section: str,
    key: str,
    unsafe_value: object,
) -> None:
    config = copy.deepcopy(canonical_config())
    config[section][key] = unsafe_value  # type: ignore[index]
    with pytest.raises(PublicationConfigurationError, match=f"{section}.{key}"):
        validate_publication_config(config)
