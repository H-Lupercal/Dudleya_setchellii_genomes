"""Validation for the approved supplementary-analysis policy."""

from __future__ import annotations

from collections.abc import Mapping


class SupplementConfigurationError(ValueError):
    """Raised when configuration departs from an approved decision-plan revision."""


def _get(config: Mapping[str, object], section: str, key: str) -> object:
    values = config.get(section)
    if not isinstance(values, Mapping) or key not in values:
        raise SupplementConfigurationError(f"Missing configuration: {section}.{key}")
    return values[key]


def validate_config(config: Mapping[str, object]) -> None:
    decision_version = _get(config, "workflow", "decision_plan_version")
    if decision_version not in {"2.5", "2.6"}:
        raise SupplementConfigurationError(f"Unsupported workflow.decision_plan_version={decision_version!r}")
    exact = {
        ("workflow", "kind"): "supplementary",
        ("workflow", "base_run_id"): "publication-20260817",
        ("likelihood_mapping", "quartets"): 100000,
        ("likelihood_mapping", "center_limit"): 0.15,
        ("likelihood_mapping", "side_trigger"): 0.20,
        ("likelihood_mapping", "split_trigger"): 0.20,
        ("seeds", "cp_tree"): 271828,
        ("seeds", "mt_tree"): 314159,
        ("seeds", "site_resampling"): 424200,
        ("seeds", "pi_resampling"): 424201,
        ("seeds", "protest"): [424210, 424211, 424212, 424213, 424214, 424215],
        ("seeds", "technical_confounders_start"): 424300,
    }
    if decision_version == "2.6":
        exact.update(
            {
                ("workflow", "run_id"): "supplement-20260824-v26",
                ("supersedes", "run_id"): "supplement-20260824",
                ("supersedes", "git_commit"): "dfbf23dacea8edd220cab88d090692e2cf7a5099",
                ("supersedes", "acceptance_sha256"): ("0682b1f8f734bff7c36d3f616628694891c9811ceecedb35d17cd665b0f51303"),
                ("supersedes", "artifact_manifest_sha256"): ("cf38c04a78ad35bf4853d0bdb94a15b64eb4b0b1bffb597e146227f5c0db90ee"),
                ("seeds", "complete_pca_protest"): [424318, 424319],
                ("seeds", "complete_pca_confounders_start"): 424320,
                ("likelihood_mapping", "mitochondria_mask_length"): 43182,
            }
        )
    scenarios = config.get("scenarios")
    if not isinstance(scenarios, Mapping):
        raise SupplementConfigurationError("Missing configuration: scenarios")
    expected_scenarios = {
        "canonical": {"dp": 5, "gq": 20, "missing": 0.20, "breadth": 0.80, "eligibility_dp": 5},
        "permissive": {"dp": 3, "gq": 15, "missing": 0.30, "breadth": 0.70, "eligibility_dp": 3},
        "strict": {"dp": 10, "gq": 30, "missing": 0.10, "breadth": 0.90, "eligibility_dp": 10},
        "mtmask70": {"mask_support": 0.70},
        "mtmask90": {"mask_support": 0.90},
    }
    errors: list[str] = []
    for (section, key), expected in exact.items():
        observed = _get(config, section, key)
        if observed != expected:
            errors.append(f"{section}.{key}={observed!r}; required {expected!r}")
    for name, values in expected_scenarios.items():
        observed = scenarios.get(name)
        if not isinstance(observed, Mapping):
            errors.append(f"scenarios.{name} missing")
            continue
        for key, expected in values.items():
            if observed.get(key) != expected:
                errors.append(f"{name}.{key}={observed.get(key)!r}; required {expected!r}")
    figures = _get(config, "figures", "families")
    if not isinstance(figures, list) or len(figures) != 6 or len(set(figures)) != 6:
        errors.append("figures.families must contain exactly six figure families")
    if errors:
        raise SupplementConfigurationError("; ".join(errors))
