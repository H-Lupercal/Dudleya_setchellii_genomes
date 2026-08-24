"""Supplementary metadata policies without modifying canonical metadata."""

from __future__ import annotations

from collections.abc import Iterable, Mapping

CORRECTION_DATE = "2026-08-24"


def apply_metadata_policy(rows: Iterable[Mapping[str, str]]) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    corrected: list[dict[str, str]] = []
    changes: list[dict[str, str]] = []
    for source in rows:
        row = dict(source)
        row["sample_analysis_eligible"] = row.get("analysis_eligible", "yes")
        row["population_inference_eligible"] = "yes"
        row["population_exclusion_reason"] = ""
        if row.get("popcode") == "DUSE":
            row["population_inference_eligible"] = "no"
            row["population_exclusion_reason"] = "population label unresolved after repository source-record audit"
            changes.append(
                {
                    "sample_id": row["sample_id"],
                    "old_popcode": "DUSE",
                    "new_popcode_or_EXCLUDED": "EXCLUDED",
                    "evidence_source": "repository source-record audit",
                    "decision_author": "supplementary pipeline policy",
                    "decision_date": CORRECTION_DATE,
                    "confidence_or_unresolved": "unresolved",
                }
            )
        corrected.append(row)
    return corrected, changes


def derive_populations(rows: Iterable[Mapping[str, str]]) -> list[dict[str, str]]:
    populations: dict[str, dict[str, str]] = {}
    for row in rows:
        if row.get("population_inference_eligible") != "yes":
            continue
        popcode = row["popcode"]
        value = {"popcode": popcode, "species": row["species"], "population_name": row["population_name"]}
        previous = populations.get(popcode)
        if previous is not None and previous != value:
            raise ValueError(f"Conflicting metadata for population {popcode}")
        populations[popcode] = value
    return [populations[key] for key in sorted(populations)]


def verification_rows() -> list[dict[str, str]]:
    return [
        {
            "entity": "DUSE",
            "issue": "population label",
            "status": "unresolved",
            "action": "exclude from population inference; retain sample-level analyses",
        },
        {
            "entity": "CY_CAS",
            "issue": "blank source species",
            "status": "supported_inference_not_independently_verified",
            "action": "retain population label; describe taxonomy inference",
        },
        {
            "entity": "CY_SIE",
            "issue": "statistical outlier",
            "status": "not_a_metadata_defect",
            "action": "retain and report as biological/statistical outlier",
        },
        {
            "entity": "TUL2",
            "issue": "conflicting display names",
            "status": "unresolved",
            "action": "do not display unsupported resolved name",
        },
    ]
