#!/usr/bin/env python3
"""Build derived figures and notes for review-response tree/PCA requests."""

from __future__ import annotations

import csv
import os
from collections import Counter
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(__file__).resolve().parents[1] / "review_response" / ".matplotlib"),
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Bio import Phylo


ROOT = Path(__file__).resolve().parents[1]
FULL_RUN = ROOT / "full_pipeline_run"
RESULTS = FULL_RUN / "results"
TREE_DIR = RESULTS / "19_bootstrap_phylogenetic_tree"
PCA_DIR = RESULTS / "15_pca"
METADATA_PATH = RESULTS / "07_downstream_sample_set" / "included_samples.tsv"
OUTDIR = ROOT / "review_response"
STRONG_SUPPORT = 95.0

GROUP_COLORS = {
    "ABAB": "#7b3294",
    "ABBE": "#c2a5cf",
    "ABMU": "#008837",
    "DUSE": "#1f78b4",
    "DUCY": "#e66101",
    "Other / legacy IDs": "#8c8c8c",
}
GROUP_ORDER = ["ABAB", "ABBE", "ABMU", "DUSE", "DUCY", "Other / legacy IDs"]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def sample_group(row: dict[str, str]) -> str:
    sample_id = row["sample_id"]
    species = row.get("species", "")
    popcode = row.get("popcode", "")
    if sample_id.startswith("ABAB_") or popcode.startswith("ABAB"):
        return "ABAB"
    if sample_id.startswith("ABBE_") or popcode.startswith("ABBE"):
        return "ABBE"
    if sample_id.startswith("ABMU_") or popcode.startswith("ABMU"):
        return "ABMU"
    if "setchellii" in species:
        return "DUSE"
    if "cymosa" in species:
        return "DUCY"
    return "Other / legacy IDs"


def load_metadata() -> dict[str, dict[str, str]]:
    metadata = {}
    for row in read_tsv(METADATA_PATH):
        row["display_group"] = sample_group(row)
        metadata[row["sample_id"]] = row
    return metadata


def outgroup_names(tree, prefixes: tuple[str, ...]) -> list[str]:
    names = []
    for terminal in tree.get_terminals():
        if terminal.name and terminal.name.startswith(prefixes):
            names.append(terminal.name)
    return sorted(names)


def root_tree(tree_path: Path, outgroup_prefixes: tuple[str, ...]):
    tree = Phylo.read(tree_path, "newick")
    names = outgroup_names(tree, outgroup_prefixes)
    if not names:
        raise RuntimeError(f"No outgroup tips matching {outgroup_prefixes} in {tree_path}")
    targets = [next(t for t in tree.get_terminals() if t.name == name) for name in names]
    tree.root_with_outgroup(targets[0], *targets[1:])
    return tree, names


def write_rooted_tree_outputs(metadata: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    rows = []
    root_specs = {
        "ABAB_ABMU": ("ABAB_", "ABMU_"),
        "ABAB": ("ABAB_",),
        "ABMU": ("ABMU_",),
    }
    for organelle in ("cpDNA", "mtDNA"):
        tree_path = TREE_DIR / f"{organelle}.primary.iqtree_ml.treefile"
        for root_label, prefixes in root_specs.items():
            tree, names = root_tree(tree_path, prefixes)
            prefix = OUTDIR / f"{organelle}.primary.rooted_{root_label}.iqtree_ml"
            newick_path = Path(f"{prefix}.treefile")
            png_path = Path(f"{prefix}.png")
            pdf_path = Path(f"{prefix}.pdf")
            svg_path = Path(f"{prefix}.svg")
            Phylo.write(tree, newick_path, "newick")
            render_tree(tree, png_path, pdf_path, svg_path, metadata, organelle, root_label)
            rows.append(
                {
                    "organelle": organelle,
                    "rooting": root_label,
                    "outgroup_tip_count": str(len(names)),
                    "outgroup_tips": ";".join(names),
                    "rooted_treefile": newick_path.as_posix(),
                    "png_path": png_path.as_posix(),
                    "pdf_path": pdf_path.as_posix(),
                    "svg_path": svg_path.as_posix(),
                }
            )
    write_tsv(
        OUTDIR / "rooted_tree_summary.tsv",
        rows,
        [
            "organelle",
            "rooting",
            "outgroup_tip_count",
            "outgroup_tips",
            "rooted_treefile",
            "png_path",
            "pdf_path",
            "svg_path",
        ],
    )
    return rows


def render_tree(tree, png_path: Path, pdf_path: Path, svg_path: Path, metadata, organelle: str, root_label: str) -> None:
    fig = plt.figure(figsize=(6, 19), constrained_layout=True)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title(
        f"{organelle} ML tree rooted with {root_label.replace('_', ' + ')}",
        fontsize=8,
        fontweight="bold",
    )
    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        label_func=lambda clade: clade.name if clade.is_terminal() else None,
        branch_labels=lambda clade: (
            str(int(round(float(clade.confidence))))
            if clade.confidence is not None
            else None
        ),
    )
    for text in ax.texts:
        label = text.get_text()
        text.set_fontsize(6)
        if label in metadata:
            text.set_color(GROUP_COLORS[metadata[label]["display_group"]])
        else:
            text.set_color("#555555")
    ax.set_xlabel("Substitutions per site", fontsize=7)
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelsize=6)
    ax.tick_params(axis="y", left=False, labelleft=False)
    ax.grid(axis="x", alpha=0.20, linewidth=0.4)
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            label=group,
            markerfacecolor=GROUP_COLORS[group],
            markersize=5,
        )
        for group in GROUP_ORDER
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=6, title="Tip group", title_fontsize=6)
    fig.savefig(png_path, dpi=400)
    fig.savefig(pdf_path)
    fig.savefig(svg_path)
    plt.close(fig)


def terminal_set(tree) -> frozenset[str]:
    return frozenset(t.name for t in tree.get_terminals() if t.name)


def canonical_split(subset: frozenset[str], all_names: frozenset[str]) -> frozenset[str]:
    complement = all_names - subset
    if len(complement) < len(subset):
        return frozenset(complement)
    if len(complement) == len(subset) and tuple(sorted(complement)) < tuple(sorted(subset)):
        return frozenset(complement)
    return subset


def strong_splits(tree, all_names: frozenset[str]) -> dict[frozenset[str], float]:
    splits = {}
    for clade in tree.get_nonterminals():
        support = clade.confidence
        if support is None or float(support) < STRONG_SUPPORT:
            continue
        subset = frozenset(t.name for t in clade.get_terminals() if t.name)
        if len(subset) < 2 or len(subset) == len(all_names):
            continue
        split = canonical_split(subset, all_names)
        if len(split) < 2:
            continue
        splits[split] = max(float(support), splits.get(split, 0.0))
    return splits


def compatible(a: frozenset[str], b: frozenset[str], all_names: frozenset[str]) -> bool:
    acomp = all_names - a
    bcomp = all_names - b
    return not (a & b and a & bcomp and acomp & b and acomp & bcomp)


def describe_split(split: frozenset[str], metadata: dict[str, dict[str, str]]) -> str:
    counts = Counter(metadata.get(name, {}).get("display_group", "Other / legacy IDs") for name in split)
    group_text = ", ".join(f"{group}={counts[group]}" for group in GROUP_ORDER if counts[group])
    examples = ", ".join(sorted(split)[:5])
    if len(split) > 5:
        examples += ", ..."
    return f"{len(split)} tips ({group_text}); examples: {examples}"


def build_tree_comparison(metadata: dict[str, dict[str, str]]) -> None:
    cptree, _ = root_tree(TREE_DIR / "cpDNA.primary.iqtree_ml.treefile", ("ABAB_", "ABMU_"))
    mttree, _ = root_tree(TREE_DIR / "mtDNA.primary.iqtree_ml.treefile", ("ABAB_", "ABMU_"))
    all_names = terminal_set(cptree) & terminal_set(mttree)
    cp = strong_splits(cptree, all_names)
    mt = strong_splits(mttree, all_names)

    shared = []
    for split in sorted(set(cp) & set(mt), key=lambda s: (-min(cp[s], mt[s]), -len(s), sorted(s))):
        shared.append(
            {
                "split_description": describe_split(split, metadata),
                "cpDNA_support": f"{cp[split]:.0f}",
                "mtDNA_support": f"{mt[split]:.0f}",
            }
        )

    conflicts = []
    seen_conflict_pairs = set()
    for source_name, source, other_name, other in (
        ("cpDNA", cp, "mtDNA", mt),
        ("mtDNA", mt, "cpDNA", cp),
    ):
        for split, support in source.items():
            if split in other:
                continue
            conflict_candidates = [
                (other_split, other_support)
                for other_split, other_support in other.items()
                if not compatible(split, other_split, all_names)
            ]
            if not conflict_candidates:
                continue
            other_split, other_support = max(conflict_candidates, key=lambda item: item[1])
            pair_key = tuple(sorted((tuple(sorted(split)), tuple(sorted(other_split)))))
            if pair_key in seen_conflict_pairs:
                continue
            seen_conflict_pairs.add(pair_key)
            conflicts.append(
                {
                    "source_tree": source_name,
                    "source_support": f"{support:.0f}",
                    "source_split": describe_split(split, metadata),
                    "conflicting_tree": other_name,
                    "conflicting_support": f"{other_support:.0f}",
                    "conflicting_split": describe_split(other_split, metadata),
                }
            )
    conflicts.sort(key=lambda row: (-(float(row["source_support"]) + float(row["conflicting_support"])), row["source_tree"]))

    write_tsv(
        OUTDIR / "strongly_supported_shared_splits.tsv",
        shared[:50],
        ["split_description", "cpDNA_support", "mtDNA_support"],
    )
    write_tsv(
        OUTDIR / "strongly_supported_conflicting_splits.tsv",
        conflicts[:50],
        [
            "source_tree",
            "source_support",
            "source_split",
            "conflicting_tree",
            "conflicting_support",
            "conflicting_split",
        ],
    )

    lines = [
        "# cpDNA vs mtDNA Tree Comparison",
        "",
        f"Strongly supported branches are defined here as UFBoot support >= {STRONG_SUPPORT:.0f}.",
        "The comparison uses the Stage 19 maximum-likelihood trees, rerooted on the ABAB + ABMU outgroup set.",
        "",
        "## Similarities",
        "",
    ]
    if shared:
        for row in shared[:8]:
            lines.append(
                f"- Shared strong split: {row['split_description']} "
                f"(cpDNA {row['cpDNA_support']}, mtDNA {row['mtDNA_support']})."
            )
    else:
        lines.append("- No identical strongly supported internal splits were found under this threshold.")
    lines.extend(["", "## Discrepancies", ""])
    if conflicts:
        for row in conflicts[:10]:
            lines.append(
                f"- {row['source_tree']} strongly supports {row['source_split']} "
                f"(support {row['source_support']}), but {row['conflicting_tree']} strongly supports an incompatible split: "
                f"{row['conflicting_split']} (support {row['conflicting_support']})."
            )
    else:
        lines.append("- No strongly supported incompatible splits were found under this threshold.")
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- The original Stage 19 trees already contain 1,000 ultrafast bootstrap replicates with BNNI correction.",
            "- These rooted trees are display/rooting derivatives of the existing Stage 19 treefiles; the ML search was not rerun.",
            "- The split tables provide the audit trail behind the short list above.",
        ]
    )
    (OUTDIR / "cpDNA_mtDNA_tree_comparison.md").write_text("\n".join(lines))


def pca_group_for_sample(sample_id: str, metadata: dict[str, dict[str, str]]) -> str:
    return metadata.get(sample_id, {}).get("display_group", "Other / legacy IDs")


def write_grouped_pca(metadata: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    rows = []
    for organelle in ("cpDNA", "mtDNA"):
        coord_path = PCA_DIR / f"{organelle}.primary.pca.coordinates.tsv"
        var_path = PCA_DIR / f"{organelle}.primary.pca.variance.tsv"
        coords = read_tsv(coord_path)
        variances = {
            row["component"]: float(row["explained_variance_ratio"])
            for row in read_tsv(var_path)
        }
        for row in coords:
            row["display_group"] = pca_group_for_sample(row["sample_id"], metadata)
        grouped_coord_path = OUTDIR / f"{organelle}.primary.pca.requested_groups.coordinates.tsv"
        write_tsv(grouped_coord_path, coords, list(coords[0].keys()))
        prefix = OUTDIR / f"{organelle}.primary.pca.requested_groups"
        png_path = Path(f"{prefix}.png")
        pdf_path = Path(f"{prefix}.pdf")
        svg_path = Path(f"{prefix}.svg")
        render_pca(coords, variances, organelle, png_path, pdf_path, svg_path)
        rows.append(
            {
                "organelle": organelle,
                "coordinates_path": grouped_coord_path.as_posix(),
                "png_path": png_path.as_posix(),
                "pdf_path": pdf_path.as_posix(),
                "svg_path": svg_path.as_posix(),
            }
        )
    write_tsv(
        OUTDIR / "pca_requested_group_summary.tsv",
        rows,
        ["organelle", "coordinates_path", "png_path", "pdf_path", "svg_path"],
    )
    return rows


def render_pca(rows, variances, organelle: str, png_path: Path, pdf_path: Path, svg_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 7), constrained_layout=True)
    for group in GROUP_ORDER:
        group_rows = [row for row in rows if row["display_group"] == group]
        if not group_rows:
            continue
        alpha = 0.82 if group != "Other / legacy IDs" else 0.35
        size = 48 if group != "Other / legacy IDs" else 24
        ax.scatter(
            [float(row["pc1"]) for row in group_rows],
            [float(row["pc2"]) for row in group_rows],
            s=size,
            alpha=alpha,
            label=f"{group} (n={len(group_rows)})",
            color=GROUP_COLORS[group],
            edgecolors="white",
            linewidths=0.35,
        )
    ax.axhline(0, color="#666666", linewidth=0.6, alpha=0.4)
    ax.axvline(0, color="#666666", linewidth=0.6, alpha=0.4)
    ax.grid(alpha=0.18, linewidth=0.5)
    ax.set_title(f"{organelle} PCA grouped as DUSE / DUCY / ABAB / ABBE / ABMU")
    ax.set_xlabel(f"PC1 ({variances['PC1'] * 100:.2f}% variance)")
    ax.set_ylabel(f"PC2 ({variances['PC2'] * 100:.2f}% variance)")
    ax.legend(title="Group", loc="best", frameon=True, fontsize=8, title_fontsize=9)
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    fig.savefig(svg_path)
    plt.close(fig)


def write_readme(rooted_rows, pca_rows) -> None:
    lines = [
        "# Review Response Outputs",
        "",
        "Derived outputs responding to the tree-rooting and PCA-legend requests.",
        "The Stage 19 IQ-TREE searches were not rerun; those trees already include 1,000 ultrafast bootstrap replicates with BNNI correction.",
        "",
        "## Publication Figures",
        "",
        "- [Paper-style figure set](publication_figures/README.md)",
        "- [Rooted cpDNA and mtDNA trees](publication_figures/figure_1_rooted_collapsed_trees.png)",
        "- [cpDNA/mtDNA branch comparison](publication_figures/figure_2_cpdna_mtdna_branch_comparison.png)",
        "- [Legended cpDNA and mtDNA PCA](publication_figures/figure_3_pca_requested_groups.png)",
        "",
        "## Rooted Trees",
        "",
    ]
    for row in rooted_rows:
        if row["rooting"] != "ABAB_ABMU":
            continue
        lines.append(
            f"- {row['organelle']} rooted with ABAB + ABMU: "
            f"[treefile]({Path(row['rooted_treefile']).name}), "
            f"[PNG]({Path(row['png_path']).name}), "
            f"[PDF]({Path(row['pdf_path']).name}), "
            f"[SVG]({Path(row['svg_path']).name})"
        )
    lines.extend(["", "Alternative ABAB-only and ABMU-only rooted versions are also included in `rooted_tree_summary.tsv`.", "", "## Tree Comparison", ""])
    lines.append("- [cpDNA_mtDNA_tree_comparison.md](cpDNA_mtDNA_tree_comparison.md)")
    lines.append("- [strongly_supported_shared_splits.tsv](strongly_supported_shared_splits.tsv)")
    lines.append("- [strongly_supported_conflicting_splits.tsv](strongly_supported_conflicting_splits.tsv)")
    lines.extend(["", "## PCA With Embedded Legend", ""])
    for row in pca_rows:
        lines.append(
            f"- {row['organelle']} PCA grouped by requested groups: "
            f"[PNG]({Path(row['png_path']).name}), "
            f"[PDF]({Path(row['pdf_path']).name}), "
            f"[SVG]({Path(row['svg_path']).name}), "
            f"[coordinates]({Path(row['coordinates_path']).name})"
        )
    lines.extend(
        [
            "",
            "Group mapping used here:",
            "",
            "- `ABAB`: ABAB-prefixed samples / D. abramsii ssp. abramsii",
            "- `ABBE`: ABBE-prefixed samples / D. abramsii ssp. bettinae",
            "- `ABMU`: ABMU-prefixed samples / D. abramsii ssp. murina",
            "- `DUSE`: samples annotated as D. setchellii",
            "- `DUCY`: samples annotated as D. cymosa",
            "- `Other / legacy IDs`: included samples without one of those labels",
            "",
        ]
    )
    (OUTDIR / "README.md").write_text("\n".join(lines))


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    metadata = load_metadata()
    rooted_rows = write_rooted_tree_outputs(metadata)
    build_tree_comparison(metadata)
    pca_rows = write_grouped_pca(metadata)
    write_readme(rooted_rows, pca_rows)
    print(f"Wrote review-response outputs to {OUTDIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
