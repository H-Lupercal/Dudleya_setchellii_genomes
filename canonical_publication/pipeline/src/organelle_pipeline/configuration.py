"""Validation for the locked publication-remediation policy."""

from __future__ import annotations

from collections.abc import Mapping


class PublicationConfigurationError(ValueError):
    """Raised when configuration departs from an acceptance-critical policy."""


def _value(config: Mapping[str, object], section: str, key: str) -> object:
    values = config.get(section)
    if not isinstance(values, Mapping) or key not in values:
        raise PublicationConfigurationError(f"Missing publication configuration: {section}.{key}")
    return values[key]


def validate_publication_config(config: Mapping[str, object]) -> None:
    """Reject settings that would no longer implement the approved remediation."""

    exact = {
        ("execution", "mapping_jobs"): 2,
        ("execution", "mapping_threads_per_job"): 8,
        ("execution", "admixture_jobs"): 4,
        ("execution", "admixture_threads_per_job"): 2,
        ("preprocessing", "qualified_quality_phred"): 20,
        ("preprocessing", "maximum_unqualified_base_percent"): 40,
        ("preprocessing", "minimum_length"): 50,
        ("preprocessing", "detect_adapters_for_pe"): True,
        ("mapping", "minimum_mapping_quality"): 20,
        ("mapping", "minimum_base_quality"): 20,
        ("mapping", "exclude_sam_flags"): 3844,
        ("mapping", "mark_duplicates"): True,
        ("mapping", "nuclear_decoy"): False,
        ("qc", "breadth_depths"): [1, 3, 5, 10],
        ("qc", "eligibility_depth"): 5,
        ("qc", "minimum_breadth"): 0.80,
        ("qc", "sample_sets"): "organelle_specific",
        ("qc", "breadth_denominator"): "organelle_unique_mappability_mask",
        ("variants", "ploidy"): 1,
        ("variants", "minimum_depth"): 5,
        ("variants", "minimum_genotype_quality"): 20,
        ("variants", "minimum_site_quality"): 30,
        ("variants", "maximum_missing_fraction"): 0.20,
        ("variants", "primary_minimum_minor_allele_count"): 1,
        ("variants", "ordination_minimum_minor_allele_count"): 2,
        ("population_genetics", "fst_estimator"): "hudson_ratio_of_sums",
        ("population_genetics", "bootstrap_block_size"): 1000,
        ("population_genetics", "bootstrap_replicates"): 1000,
        ("phylogeny", "model"): "MFP",
        ("phylogeny", "shalrt_replicates"): 1000,
        ("phylogeny", "ultrafast_bootstrap_replicates"): 1000,
        ("phylogeny", "bnni"): True,
        ("phylogeny", "rooting"): "unrooted",
        ("admixture", "role"): "supplementary",
        ("admixture", "minimum_k"): 1,
        ("admixture", "maximum_k"): 12,
        ("admixture", "replicates"): 10,
        ("concatenation", "role"): "supplementary",
        ("concatenation", "sample_set"): "shared_intersection",
        ("concatenation", "partitioned"): True,
    }
    mismatches = [
        f"{section}.{key}={_value(config, section, key)!r} (required {expected!r})"
        for (section, key), expected in exact.items()
        if _value(config, section, key) != expected
    ]
    required_flags = 4 | 256 | 512 | 1024 | 2048
    observed_flags = int(_value(config, "mapping", "exclude_sam_flags"))
    if observed_flags & required_flags != required_flags:
        mismatches.append(f"mapping.exclude_sam_flags={observed_flags!r} does not include {required_flags}")
    if mismatches:
        raise PublicationConfigurationError("Configuration departs from the approved publication policy: " + "; ".join(mismatches))
