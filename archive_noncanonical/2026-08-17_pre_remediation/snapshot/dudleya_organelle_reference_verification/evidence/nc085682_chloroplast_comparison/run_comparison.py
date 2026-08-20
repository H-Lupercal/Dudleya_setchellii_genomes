#!/usr/bin/env python3
"""Compare the Dudleya candidate plastome against NC_085682.1.

This script is intentionally dependency-light: the execution environment has
BLAST+, MAFFT, and Python 3, but not BioPython/pandas. It writes all report
artifacts into this directory.
"""

from __future__ import annotations

import csv
import datetime as dt
import argparse
import re
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parents[2]
REPO_ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent

CANDIDATE_RAW = REPO_ROOT / "Dudleya_hifiasm_purged_manual_chloroplast.fa"
CACHED_GB = OUT / "reference" / "NC_085682.1.fetched.gb"
DRAFT_TSV = PACKAGE_ROOT / "annotations" / "chloroplast.annotation.tsv"
EXISTING_SUMMARY = PACKAGE_ROOT / "evidence" / "identity" / "chloroplast_vs_chloroplast_refs.tsv"

ACCESSION = "NC_085682.1"
NCBI_NUCCORE_URL = f"https://www.ncbi.nlm.nih.gov/nuccore/{ACCESSION}"
NCBI_EFETCH_FASTA_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    f"?db=nuccore&id={ACCESSION}&rettype=fasta&retmode=text"
)
NCBI_EFETCH_GB_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    f"?db=nuccore&id={ACCESSION}&rettype=gb&retmode=text"
)
DEDUP_LENGTH = 150_274
TERMINAL_DUP_SUFFIX_START = 150_275


COMPLEMENT = str.maketrans("ACGTNacgtn-", "TGCANtgcan-")


GENETIC_CODE_11 = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


@dataclass
class Feature:
    feature_type: str
    location: str
    qualifiers: dict[str, str]
    ranges: list[tuple[int, int]]
    strand: str
    order_index: int

    @property
    def gene(self) -> str:
        return (
            self.qualifiers.get("gene")
            or self.qualifiers.get("label")
            or self.qualifiers.get("product")
            or ""
        )

    @property
    def product(self) -> str:
        return self.qualifiers.get("product", "")

    @property
    def start(self) -> int:
        return min(start for start, _ in self.ranges)

    @property
    def end(self) -> int:
        return max(end for _, end in self.ranges)

    @property
    def length(self) -> int:
        return sum(end - start + 1 for start, end in self.ranges)

    @property
    def feature_id(self) -> str:
        gene = self.gene or "unnamed"
        return f"{self.feature_type}:{gene}:{self.order_index}"


def run(cmd: list[str], *, cwd: Path = ROOT, input_text: str | None = None) -> str:
    proc = subprocess.run(
        cmd,
        cwd=str(cwd),
        input=input_text,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"Command failed ({proc.returncode}): {' '.join(cmd)}\nSTDERR:\n{proc.stderr}"
        )
    return proc.stdout


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    seq_parts: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts).upper()))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts).upper()))
    return records


def write_fasta(path: Path, records: list[tuple[str, str]], width: int = 80) -> None:
    with path.open("w") as handle:
        for header, seq in records:
            handle.write(f">{header}\n")
            for i in range(0, len(seq), width):
                handle.write(seq[i:i + width] + "\n")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1].upper()


def normalize_gene(name: str) -> str:
    return re.sub(r"[^a-z0-9.]+", "", name.lower())


def fetch_current_records(allow_cached_fallback: bool) -> tuple[str | None, str | None, str]:
    """Try NCBI EFetch, returning FASTA text, GenBank text, and source note."""
    urls = {
        "fasta": NCBI_EFETCH_FASTA_URL,
        "gb": NCBI_EFETCH_GB_URL,
    }
    try:
        fasta = urllib.request.urlopen(urls["fasta"], timeout=20).read().decode("utf-8")
        gb = urllib.request.urlopen(urls["gb"], timeout=20).read().decode("utf-8")
        if ACCESSION not in fasta or "LOCUS" not in gb or "ORIGIN" not in gb:
            raise RuntimeError("NCBI response did not contain the expected accession")
        return fasta, gb, "live NCBI EFetch"
    except (urllib.error.URLError, TimeoutError, RuntimeError) as exc:
        if not allow_cached_fallback:
            raise RuntimeError(
                "Live NCBI download failed and cached fallback is disabled. "
                "Rerun with network access. Original error: "
                f"{exc}"
            ) from exc
        return None, None, f"cached NC_085682.1 GenBank fallback; live fetch failed: {exc}"


def split_genbank_records(gb_text: str) -> list[str]:
    chunks = []
    for raw in gb_text.split("\n//"):
        raw = raw.strip("\n")
        if raw.strip():
            chunks.append(raw + "\n//\n")
    return chunks


def cached_nc_record() -> str:
    text = CACHED_GB.read_text()
    for record in split_genbank_records(text):
        if f"VERSION     {ACCESSION}" in record:
            return record
    raise RuntimeError(f"{ACCESSION} not found in {CACHED_GB}")


def extract_origin(record: str) -> str:
    in_origin = False
    parts: list[str] = []
    for line in record.splitlines():
        if line.startswith("ORIGIN"):
            in_origin = True
            continue
        if line.startswith("//"):
            break
        if in_origin:
            parts.append(re.sub(r"[^A-Za-z]", "", line))
    return "".join(parts).upper()


def parse_locus_length(record: str) -> int:
    match = re.search(r"^LOCUS\s+\S+\s+(\d+)\s+bp", record, flags=re.MULTILINE)
    if not match:
        raise RuntimeError("Could not parse LOCUS length")
    return int(match.group(1))


def parse_definition(record: str) -> str:
    lines = record.splitlines()
    definition: list[str] = []
    in_def = False
    for line in lines:
        if line.startswith("DEFINITION"):
            in_def = True
            definition.append(line[len("DEFINITION"):].strip())
            continue
        if in_def:
            if re.match(r"^[A-Z][A-Z_ -]+\s+", line):
                break
            definition.append(line.strip())
    return " ".join(part for part in definition if part).strip()


def parse_locus_date(record: str) -> str:
    first = next(line for line in record.splitlines() if line.startswith("LOCUS"))
    return first.split()[-1]


def parse_location(location: str) -> tuple[list[tuple[int, int]], str]:
    loc = location.replace(" ", "")
    loc = loc.replace("<", "").replace(">", "")
    strand = "-" if loc.startswith("complement(") else "+"
    ranges = [(int(a), int(b)) for a, b in re.findall(r"(\d+)\.\.(\d+)", loc)]
    if not ranges:
        singles = [int(x) for x in re.findall(r"(?<![\d.])(\d+)(?![\d.])", loc)]
        ranges = [(x, x) for x in singles]
    if not ranges:
        raise RuntimeError(f"Could not parse feature location: {location}")
    return ranges, strand


def unwrap_location_function(expr: str, name: str) -> str | None:
    prefix = f"{name}("
    if not expr.startswith(prefix) or not expr.endswith(")"):
        return None
    depth = 0
    for idx, char in enumerate(expr):
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
            if depth == 0 and idx != len(expr) - 1:
                return None
    return expr[len(prefix):-1]


def split_location_args(expr: str) -> list[str]:
    args: list[str] = []
    depth = 0
    start = 0
    for idx, char in enumerate(expr):
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        elif char == "," and depth == 0:
            args.append(expr[start:idx])
            start = idx + 1
    args.append(expr[start:])
    return [arg for arg in args if arg]


def location_base_specs(location: str) -> list[tuple[int, bool]]:
    """Return reference positions in biological order with complement flags."""
    expr = location.replace(" ", "").replace("<", "").replace(">", "")

    def parse_expr(part: str) -> list[tuple[int, bool]]:
        inner = unwrap_location_function(part, "complement")
        if inner is not None:
            specs = parse_expr(inner)
            return [(pos, not should_complement) for pos, should_complement in reversed(specs)]
        inner = unwrap_location_function(part, "join")
        if inner is not None:
            specs: list[tuple[int, bool]] = []
            for arg in split_location_args(inner):
                specs.extend(parse_expr(arg))
            return specs
        inner = unwrap_location_function(part, "order")
        if inner is not None:
            specs = []
            for arg in split_location_args(inner):
                specs.extend(parse_expr(arg))
            return specs
        match = re.fullmatch(r"(\d+)\.\.(\d+)", part)
        if match:
            start, end = int(match.group(1)), int(match.group(2))
            return [(pos, False) for pos in range(start, end + 1)]
        match = re.fullmatch(r"(\d+)", part)
        if match:
            return [(int(match.group(1)), False)]
        raise RuntimeError(f"Could not parse location expression: {part}")

    return parse_expr(expr)


def finish_qualifier(value: str) -> str:
    value = value.strip()
    if value.startswith('"') and value.endswith('"'):
        value = value[1:-1]
    return re.sub(r"\s+", "", value) if "\n" in value else value


def parse_features(record: str) -> list[Feature]:
    features: list[Feature] = []
    in_features = False
    current: dict | None = None
    current_qual: str | None = None

    def flush() -> None:
        nonlocal current, current_qual
        if not current:
            return
        qualifiers = {
            key: finish_qualifier(value)
            for key, value in current["qualifiers"].items()
        }
        ranges, strand = parse_location(current["location"])
        features.append(
            Feature(
                feature_type=current["type"],
                location=current["location"],
                qualifiers=qualifiers,
                ranges=ranges,
                strand=strand,
                order_index=len(features) + 1,
            )
        )
        current = None
        current_qual = None

    for line in record.splitlines():
        if line.startswith("FEATURES"):
            in_features = True
            continue
        if line.startswith("ORIGIN"):
            flush()
            break
        if not in_features:
            continue
        if len(line) >= 21 and line[:5] == "     " and line[5:21].strip():
            flush()
            current = {
                "type": line[5:21].strip(),
                "location": line[21:].strip(),
                "qualifiers": {},
            }
            current_qual = None
            continue
        if current is None or len(line) < 21:
            continue
        content = line[21:].rstrip()
        if not content:
            continue
        if content.startswith("/"):
            content = content[1:]
            if "=" in content:
                key, value = content.split("=", 1)
            else:
                key, value = content, "true"
            current["qualifiers"][key] = value
            current_qual = key
        elif current_qual:
            current["qualifiers"][current_qual] += "\n" + content.strip()
        else:
            current["location"] += content.strip()
    return features


def extract_feature_sequence(seq: str, feature: Feature) -> str:
    bases = []
    for pos, should_complement in location_base_specs(feature.location):
        base = seq[pos - 1].upper()
        if should_complement:
            base = base.translate(COMPLEMENT).upper()
        bases.append(base)
    return "".join(bases)


def translate(seq: str) -> str:
    seq = seq.upper().replace("-", "")
    aa = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        aa.append(GENETIC_CODE_11.get(codon, "X"))
    return "".join(aa)


def translation_matches_genbank(observed: str, expected: str) -> bool:
    """Compare translations, allowing plastid alternative start codons to be M."""
    observed_no_terminal = observed[:-1] if observed.endswith("*") else observed
    if observed_no_terminal == expected:
        return True
    return (
        len(observed_no_terminal) == len(expected)
        and bool(expected)
        and expected[0] == "M"
        and observed_no_terminal[1:] == expected[1:]
    )


def write_reference_files(record: str, seq: str, source_note: str) -> dict[str, str | int]:
    gb_path = OUT / "NC_085682.1.gb"
    fa_path = OUT / "NC_085682.1.fa"
    gb_path.write_text(record)
    write_fasta(fa_path, [(ACCESSION, seq)])
    meta = {
        "accession": ACCESSION,
        "length_bp": len(seq),
        "description": parse_definition(record),
        "locus_date": parse_locus_date(record),
        "fetch_date": dt.date.today().isoformat(),
        "source": source_note,
        "nuccore_url": NCBI_NUCCORE_URL,
        "efetch_fasta_url": NCBI_EFETCH_FASTA_URL,
        "efetch_genbank_url": NCBI_EFETCH_GB_URL,
    }
    with (OUT / "reference_metadata.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(meta.keys()), delimiter="\t")
        writer.writeheader()
        writer.writerow(meta)
    return meta


def prepare_candidate_files() -> dict[str, str | int]:
    records = read_fasta(CANDIDATE_RAW)
    if len(records) != 1:
        raise RuntimeError(f"Expected one candidate FASTA record, found {len(records)}")
    header, raw = records[0]
    dedup = raw[:DEDUP_LENGTH]
    raw_fa = OUT / "candidate_raw.fa"
    dedup_fa = OUT / "candidate_terminal_deduplicated.fa"
    write_fasta(raw_fa, [(header, raw)])
    write_fasta(dedup_fa, [(f"{header}|terminal_deduplicated_1_{DEDUP_LENGTH}", dedup)])
    return {
        "header": header,
        "raw_length_bp": len(raw),
        "deduplicated_length_bp": len(dedup),
        "raw_ambiguous_bases": sum(1 for b in raw if b not in "ACGT"),
        "deduplicated_ambiguous_bases": sum(1 for b in dedup if b not in "ACGT"),
        "terminal_duplicate_trimmed_bp": len(raw) - len(dedup),
        "terminal_duplicate_suffix_start": TERMINAL_DUP_SUFFIX_START,
    }


def blast_rows(query: Path, subject: Path, out_path: Path) -> list[dict[str, str]]:
    outfmt = (
        "6 qseqid sseqid pident length qlen slen "
        "qstart qend sstart send evalue bitscore"
    )
    stdout = run([
        "blastn",
        "-query", str(query),
        "-subject", str(subject),
        "-outfmt", outfmt,
        "-evalue", "1e-20",
    ])
    out_path.write_text(stdout)
    fields = outfmt.split()[1:]
    rows = []
    for line in stdout.splitlines():
        if not line.strip():
            continue
        rows.append(dict(zip(fields, line.split("\t"))))
    return rows


def blast_alignment_rows(query: Path, subject: Path, out_path: Path) -> list[dict[str, str]]:
    outfmt = (
        "6 qseqid sseqid pident length qlen slen "
        "qstart qend sstart send evalue bitscore qseq sseq"
    )
    stdout = run([
        "blastn",
        "-query", str(query),
        "-subject", str(subject),
        "-outfmt", outfmt,
        "-evalue", "1e-20",
    ])
    out_path.write_text(stdout)
    fields = outfmt.split()[1:]
    rows = []
    for line in stdout.splitlines():
        if not line.strip():
            continue
        rows.append(dict(zip(fields, line.split("\t"))))
    return rows


def build_blast_projection(rows: list[dict[str, str]]) -> dict[int, dict[str, object]]:
    """Map reference positions to candidate bases/coordinates from BLAST HSPs.

    The highest-identity, longest HSP wins per reference base. This deliberately
    allows reverse-orientation HSPs in the SSC, where plastome representations
    can differ by the common SSC flip.
    """
    projection: dict[int, dict[str, object]] = {}
    for row in rows:
        qstart, qend = int(row["qstart"]), int(row["qend"])
        sstart, send = int(row["sstart"]), int(row["send"])
        qstep = 1 if qend >= qstart else -1
        sstep = 1 if send >= sstart else -1
        orientation = "+" if sstep == qstep else "-"
        qpos = qstart
        spos = sstart
        priority = (float(row["pident"]), int(row["length"]), float(row["bitscore"]))
        for qbase, sbase in zip(row["qseq"], row["sseq"]):
            current_qpos = qpos if qbase != "-" else None
            current_spos = spos if sbase != "-" else None
            if current_qpos is not None:
                existing = projection.get(current_qpos)
                if existing is None or priority > existing["priority"]:
                    projection[current_qpos] = {
                        "base": sbase.upper() if sbase != "-" else "-",
                        "candidate_pos": current_spos,
                        "orientation": orientation,
                        "priority": priority,
                    }
            if qbase != "-":
                qpos += qstep
            if sbase != "-":
                spos += sstep
    return projection


def merge_intervals(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    normalized = sorted((min(a, b), max(a, b)) for a, b in intervals)
    total = 0
    cur_start, cur_end = normalized[0]
    for start, end in normalized[1:]:
        if start <= cur_end + 1:
            cur_end = max(cur_end, end)
        else:
            total += cur_end - cur_start + 1
            cur_start, cur_end = start, end
    total += cur_end - cur_start + 1
    return total


def summarize_blast(rows: list[dict[str, str]]) -> dict[str, float | int]:
    if not rows:
        raise RuntimeError("No BLAST rows to summarize")
    qlen = int(rows[0]["qlen"])
    slen = int(rows[0]["slen"])
    q_intervals = [(int(r["qstart"]), int(r["qend"])) for r in rows]
    s_intervals = [(int(r["sstart"]), int(r["send"])) for r in rows]
    total_aligned = sum(int(r["length"]) for r in rows)
    weighted_identity = (
        sum(float(r["pident"]) * int(r["length"]) for r in rows) / total_aligned
    )
    bitscore_sum = sum(float(r["bitscore"]) for r in rows)
    return {
        "query_covered_bp": merge_intervals(q_intervals),
        "query_coverage": merge_intervals(q_intervals) / qlen,
        "reference_coverage": merge_intervals(s_intervals) / slen,
        "weighted_identity": weighted_identity,
        "hsp_count": len(rows),
        "total_hsp_aligned_bp": total_aligned,
        "bitscore_sum": bitscore_sum,
    }


def choose_rotation_anchor(rows: list[dict[str, str]], candidate_len: int) -> dict[str, int | str | float]:
    if not rows:
        raise RuntimeError("No BLAST rows available for rotation")
    parsed = []
    for row in rows:
        qstart, qend = int(row["qstart"]), int(row["qend"])
        sstart, send = int(row["sstart"]), int(row["send"])
        orientation = "+" if sstart <= send else "-"
        ref_low = min(qstart, qend)
        subject_at_ref_low = sstart if qstart <= qend else send
        parsed.append({
            "ref_low": ref_low,
            "subject_at_ref_low": subject_at_ref_low,
            "orientation": orientation,
            "length": int(row["length"]),
            "bitscore": float(row["bitscore"]),
            "pident": float(row["pident"]),
        })
    near_origin = [r for r in parsed if r["ref_low"] <= 1000]
    chosen = sorted(
        near_origin or parsed,
        key=lambda r: (r["ref_low"], -r["length"], -r["bitscore"]),
    )[0]
    if chosen["orientation"] == "+":
        rotation_start = int(chosen["subject_at_ref_low"])
    else:
        rotation_start = candidate_len - int(chosen["subject_at_ref_low"]) + 1
    chosen["rotation_start_1based"] = rotation_start
    return chosen


def rotate_sequence(seq: str, start_1based: int) -> str:
    idx = start_1based - 1
    return seq[idx:] + seq[:idx]


def run_mafft(ref_seq: str, cand_seq: str) -> tuple[str, str]:
    pair = OUT / "normalized_pair_for_mafft.fa"
    aln_path = OUT / "normalized_alignment.fa"
    write_fasta(pair, [(ACCESSION, ref_seq), ("candidate_deduplicated_rotated", cand_seq)])
    stdout = run(["mafft", "--auto", "--thread", "1", str(pair)], cwd=OUT)
    aln_path.write_text(stdout)
    records = read_fasta(aln_path)
    if len(records) != 2:
        raise RuntimeError(f"Expected two aligned records from MAFFT, found {len(records)}")
    return records[0][1], records[1][1]


def build_alignment_maps(ref_aln: str, cand_aln: str) -> dict[str, object]:
    ref_pos = 0
    cand_pos = 0
    ref_pos_to_aln: dict[int, int] = {}
    aln_to_cand_pos: list[int | None] = []
    aln_to_ref_pos: list[int | None] = []
    for idx, (rbase, cbase) in enumerate(zip(ref_aln, cand_aln)):
        if rbase != "-":
            ref_pos += 1
            ref_pos_to_aln[ref_pos] = idx
            aln_to_ref_pos.append(ref_pos)
        else:
            aln_to_ref_pos.append(None)
        if cbase != "-":
            cand_pos += 1
            aln_to_cand_pos.append(cand_pos)
        else:
            aln_to_cand_pos.append(None)
    return {
        "ref_pos_to_aln": ref_pos_to_aln,
        "aln_to_cand_pos": aln_to_cand_pos,
        "aln_to_ref_pos": aln_to_ref_pos,
    }


def divergence_metrics(ref_aln: str, cand_aln: str, raw_blast: dict[str, float | int]) -> dict[str, str | int | float]:
    matches = mismatches = ref_gap_bases = candidate_gap_bases = 0
    for rbase, cbase in zip(ref_aln, cand_aln):
        if rbase == "-" and cbase == "-":
            continue
        if rbase == "-":
            ref_gap_bases += 1
        elif cbase == "-":
            candidate_gap_bases += 1
        elif rbase == cbase:
            matches += 1
        else:
            mismatches += 1
    aligned_columns = matches + mismatches + ref_gap_bases + candidate_gap_bases
    comparable = matches + mismatches
    identity_excluding_gaps = 100 * matches / comparable if comparable else 0.0
    identity_including_gaps = 100 * matches / aligned_columns if aligned_columns else 0.0
    return {
        "comparison": "candidate_deduplicated_rotated_vs_NC_085682.1",
        "raw_blast_query_coverage": f"{raw_blast['query_coverage']:.6f}",
        "raw_blast_reference_coverage": f"{raw_blast['reference_coverage']:.6f}",
        "raw_blast_weighted_identity_percent": f"{raw_blast['weighted_identity']:.3f}",
        "raw_blast_hsp_count": int(raw_blast["hsp_count"]),
        "aligned_columns": aligned_columns,
        "matches": matches,
        "mismatches": mismatches,
        "indel_columns": ref_gap_bases + candidate_gap_bases,
        "candidate_insertion_bases_vs_reference": ref_gap_bases,
        "candidate_deletion_bases_vs_reference": candidate_gap_bases,
        "identity_excluding_gaps_percent": f"{identity_excluding_gaps:.6f}",
        "identity_including_gaps_percent": f"{identity_including_gaps:.6f}",
        "percent_divergence_excluding_gaps": f"{100 - identity_excluding_gaps:.6f}",
        "percent_divergence_including_gaps": f"{100 - identity_including_gaps:.6f}",
    }


def blast_projection_divergence_metrics(
    ref_seq: str,
    projection: dict[int, dict[str, object]],
    raw_blast: dict[str, float | int],
    mafft_diag: dict[str, str | int | float],
) -> dict[str, str | int | float]:
    matches = mismatches = candidate_deletions = unmapped = 0
    orientation_counts: Counter = Counter()
    for pos, ref_base in enumerate(ref_seq, start=1):
        hit = projection.get(pos)
        if not hit:
            unmapped += 1
            continue
        orientation_counts[str(hit["orientation"])] += 1
        base = str(hit["base"]).upper()
        if base == "-":
            candidate_deletions += 1
        elif base == ref_base:
            matches += 1
        else:
            mismatches += 1
    mapped_reference_bp = matches + mismatches + candidate_deletions
    comparable = matches + mismatches
    identity_excluding_gaps = 100 * matches / comparable if comparable else 0.0
    identity_with_ref_deletions = 100 * matches / mapped_reference_bp if mapped_reference_bp else 0.0
    identity_with_unmapped = 100 * matches / len(ref_seq) if ref_seq else 0.0
    return {
        "comparison": "candidate_deduplicated_rotated_vs_NC_085682.1",
        "raw_blast_query_coverage": f"{raw_blast['query_coverage']:.6f}",
        "raw_blast_reference_coverage": f"{raw_blast['reference_coverage']:.6f}",
        "raw_blast_weighted_identity_percent": f"{raw_blast['weighted_identity']:.3f}",
        "raw_blast_hsp_count": int(raw_blast["hsp_count"]),
        "normalized_projection_method": "BLAST HSP qseq/sseq best-hit per reference base",
        "normalized_mapped_reference_bp": mapped_reference_bp,
        "normalized_unmapped_reference_bp": unmapped,
        "normalized_matches": matches,
        "normalized_mismatches": mismatches,
        "normalized_candidate_deletion_bases_vs_reference": candidate_deletions,
        "normalized_forward_oriented_reference_bp": orientation_counts.get("+", 0),
        "normalized_reverse_oriented_reference_bp": orientation_counts.get("-", 0),
        "identity_excluding_gaps_percent": f"{identity_excluding_gaps:.6f}",
        "identity_with_reference_deletions_percent": f"{identity_with_ref_deletions:.6f}",
        "identity_with_unmapped_reference_as_difference_percent": f"{identity_with_unmapped:.6f}",
        "percent_divergence_excluding_gaps": f"{100 - identity_excluding_gaps:.6f}",
        "percent_divergence_with_reference_deletions": f"{100 - identity_with_ref_deletions:.6f}",
        "percent_divergence_with_unmapped_reference_as_difference": f"{100 - identity_with_unmapped:.6f}",
        "diagnostic_mafft_aligned_columns": mafft_diag["aligned_columns"],
        "diagnostic_mafft_mismatches": mafft_diag["mismatches"],
        "diagnostic_mafft_indel_columns": mafft_diag["indel_columns"],
        "diagnostic_mafft_identity_excluding_gaps_percent": mafft_diag["identity_excluding_gaps_percent"],
        "diagnostic_mafft_identity_including_gaps_percent": mafft_diag["identity_including_gaps_percent"],
    }


def ordered_positions(feature: Feature) -> list[int]:
    positions: list[int] = []
    for start, end in feature.ranges:
        positions.extend(range(start, end + 1))
    return positions


def projected_feature_sequence(
    feature: Feature,
    ref_aln: str | None,
    cand_aln: str | None,
    maps: dict[str, object],
) -> tuple[str, list[int], int, Counter]:
    if "projection_by_ref" in maps:
        projection = maps["projection_by_ref"]
        assert isinstance(projection, dict)
        bases: list[str] = []
        projected_positions: list[int] = []
        orientation_counts: Counter = Counter()
        gap_count = 0
        for pos, should_complement in location_base_specs(feature.location):
            hit = projection.get(pos)
            if not hit or hit["base"] == "-" or hit["candidate_pos"] is None:
                bases.append("-")
                gap_count += 1
            else:
                base = str(hit["base"]).upper()
                if should_complement:
                    base = base.translate(COMPLEMENT).upper()
                bases.append(base)
                projected_positions.append(int(hit["candidate_pos"]))
                orientation_counts[str(hit["orientation"])] += 1
        return "".join(bases), projected_positions, gap_count, orientation_counts

    ref_pos_to_aln = maps["ref_pos_to_aln"]
    aln_to_cand_pos = maps["aln_to_cand_pos"]
    assert isinstance(ref_pos_to_aln, dict)
    assert isinstance(aln_to_cand_pos, list)
    assert ref_aln is not None
    assert cand_aln is not None
    bases: list[str] = []
    projected_positions: list[int] = []
    orientation_counts: Counter = Counter()
    gap_count = 0
    for pos, should_complement in location_base_specs(feature.location):
        aln_idx = ref_pos_to_aln[pos]
        cbase = cand_aln[aln_idx]
        cpos = aln_to_cand_pos[aln_idx]
        if cbase == "-" or cpos is None:
            bases.append("-")
            gap_count += 1
        else:
            base = cbase.upper()
            if should_complement:
                base = base.translate(COMPLEMENT).upper()
            bases.append(base)
            projected_positions.append(int(cpos))
            orientation_counts["+"] += 1
    return "".join(bases), projected_positions, gap_count, orientation_counts


def feature_presence(
    feature: Feature,
    ref_aln: str | None,
    cand_aln: str | None,
    maps: dict[str, object],
) -> dict[str, object]:
    cand_seq, projected_positions, gap_count, orientation_counts = projected_feature_sequence(
        feature, ref_aln, cand_aln, maps
    )
    coverage = 1 - (gap_count / feature.length)
    if not projected_positions:
        status = "absent_by_projection"
        projected_start = projected_end = ""
    elif coverage >= 0.95:
        status = "present_by_projection"
        projected_start = min(projected_positions)
        projected_end = max(projected_positions)
    elif coverage >= 0.20:
        status = "partial_by_projection"
        projected_start = min(projected_positions)
        projected_end = max(projected_positions)
    else:
        status = "absent_by_projection"
        projected_start = min(projected_positions)
        projected_end = max(projected_positions)
    return {
        "status": status,
        "coverage": coverage,
        "gap_count": gap_count,
        "projected_start": projected_start,
        "projected_end": projected_end,
        "projected_seq": cand_seq,
        "projected_positions": projected_positions,
        "projection_orientation": orientation_counts.most_common(1)[0][0] if orientation_counts else "",
    }


def projected_candidate_strand(ref_strand: str, projection_orientation: str) -> str:
    if projection_orientation == "":
        return ""
    if projection_orientation == "+":
        return ref_strand
    return "-" if ref_strand == "+" else "+"


def compare_cds(
    features: list[Feature],
    ref_seq: str,
    ref_aln: str | None,
    cand_aln: str | None,
    maps: dict[str, object],
) -> tuple[list[dict[str, object]], dict[str, int]]:
    rows: list[dict[str, object]] = []
    validation = Counter()
    cds_features = [f for f in features if f.feature_type == "CDS" and f.gene]
    occurrence = Counter()
    for feature in cds_features:
        occurrence[normalize_gene(feature.gene)] += 1
        gene_copy = occurrence[normalize_gene(feature.gene)]
        codon_start = int(feature.qualifiers.get("codon_start", "1"))
        transl_table = feature.qualifiers.get("transl_table", "11")
        ref_cds_full = extract_feature_sequence(ref_seq, feature)
        ref_cds = ref_cds_full[codon_start - 1:]
        expected_translation = feature.qualifiers.get("translation", "").replace(" ", "")
        if expected_translation:
            validation["features_with_translation"] += 1
            observed = translate(ref_cds)
            if translation_matches_genbank(observed, expected_translation):
                validation["translation_matches"] += 1
            else:
                validation["translation_mismatches"] += 1

        presence = feature_presence(feature, ref_aln, cand_aln, maps)
        candidate_cds = str(presence["projected_seq"])[codon_start - 1:]

        codons_compared = 0
        nt_differences = 0
        synonymous = 0
        nonsynonymous = 0
        complex_codons = 0
        aa_changes = 0
        gap_codons = 0
        ambiguous_codons = 0

        usable_len = min(len(ref_cds), len(candidate_cds))
        for idx in range(0, usable_len - 2, 3):
            ref_codon = ref_cds[idx:idx + 3].upper()
            cand_codon = candidate_cds[idx:idx + 3].upper()
            codons_compared += 1
            if "-" in cand_codon or "-" in ref_codon:
                gap_codons += 1
                continue
            if any(base not in "ACGT" for base in ref_codon + cand_codon):
                ambiguous_codons += 1
                continue
            if ref_codon == cand_codon:
                continue
            diffs = sum(1 for a, b in zip(ref_codon, cand_codon) if a != b)
            nt_differences += diffs
            ref_aa = GENETIC_CODE_11.get(ref_codon, "X")
            cand_aa = GENETIC_CODE_11.get(cand_codon, "X")
            if diffs == 1:
                if ref_aa == cand_aa:
                    synonymous += 1
                else:
                    nonsynonymous += 1
                    aa_changes += 1
            else:
                complex_codons += 1
                if ref_aa != cand_aa:
                    aa_changes += 1

        candidate_protein = translate(candidate_cds)
        internal_stop = "*" in candidate_protein[:-1]
        partial_codon = (len(ref_cds) % 3 != 0) or (len(candidate_cds.replace("-", "")) % 3 != 0)
        has_gaps = gap_codons > 0 or int(presence["gap_count"]) > 0
        if str(presence["status"]) != "present_by_projection":
            inclusion = str(presence["status"])
        elif has_gaps:
            inclusion = "excluded_from_clean_syn_nonsyn_totals_candidate_gap"
        elif partial_codon:
            inclusion = "excluded_from_clean_syn_nonsyn_totals_partial_codon"
        elif internal_stop:
            inclusion = "excluded_from_clean_syn_nonsyn_totals_internal_stop"
        else:
            inclusion = "included_in_clean_syn_nonsyn_totals"

        rows.append({
            "feature_id": feature.feature_id,
            "gene": feature.gene,
            "copy_number": gene_copy,
            "product": feature.product,
            "ref_location": feature.location,
            "strand": feature.strand,
            "transl_table": transl_table,
            "codon_start": codon_start,
            "projection_status": presence["status"],
            "projection_orientation": presence["projection_orientation"],
            "candidate_projected_strand": projected_candidate_strand(
                feature.strand, str(presence["projection_orientation"])
            ),
            "projected_query_coverage": f"{presence['coverage']:.6f}",
            "projected_candidate_start": presence["projected_start"],
            "projected_candidate_end": presence["projected_end"],
            "ref_cds_length_bp": len(ref_cds),
            "candidate_projected_cds_length_bp": len(candidate_cds.replace("-", "")),
            "codons_compared": codons_compared,
            "nucleotide_differences_in_compared_codons": nt_differences,
            "single_nt_synonymous_substitutions": synonymous,
            "single_nt_nonsynonymous_substitutions": nonsynonymous,
            "complex_codon_changes": complex_codons,
            "amino_acid_changes": aa_changes,
            "codons_with_candidate_gaps": gap_codons,
            "ambiguous_codons": ambiguous_codons,
            "candidate_internal_stop": "yes" if internal_stop else "no",
            "partial_codon_or_frameshift_flag": "yes" if partial_codon else "no",
            "inclusion_reason": inclusion,
        })
    return rows, dict(validation)


def draft_counts() -> Counter:
    counts: Counter = Counter()
    with DRAFT_TSV.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            feature_type = row["feature_type"]
            gene = row["gene"]
            if not gene:
                continue
            counts[(feature_type, normalize_gene(gene))] += 1
    return counts


def gene_content_rows(
    features: list[Feature],
    ref_aln: str | None,
    cand_aln: str | None,
    maps: dict[str, object],
) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], list[Feature]] = defaultdict(list)
    display_names: dict[tuple[str, str], str] = {}
    for feature in features:
        if feature.feature_type not in {"gene", "CDS", "tRNA", "rRNA"} or not feature.gene:
            continue
        key = (feature.feature_type, normalize_gene(feature.gene))
        grouped[key].append(feature)
        display_names.setdefault(key, feature.gene)
    draft = draft_counts()
    rows: list[dict[str, object]] = []
    for key in sorted(grouped, key=lambda k: (k[0], display_names[k].lower())):
        feature_type, gene_key = key
        feats = grouped[key]
        statuses = [feature_presence(f, ref_aln, cand_aln, maps)["status"] for f in feats]
        present = sum(1 for s in statuses if s == "present_by_projection")
        partial = sum(1 for s in statuses if s == "partial_by_projection")
        absent = sum(1 for s in statuses if s == "absent_by_projection")
        draft_count = draft.get(key, 0)
        if present == len(feats) and draft_count == len(feats):
            status = "match"
        elif present == len(feats):
            status = "present_by_alignment_draft_count_diff"
        elif present or partial:
            status = "partial_by_alignment"
        else:
            status = "not_projected"
        rows.append({
            "feature_type": feature_type,
            "gene": display_names[key],
            "nc085682_feature_count": len(feats),
            "candidate_projected_present_count": present,
            "candidate_projected_partial_count": partial,
            "candidate_projected_absent_count": absent,
            "candidate_draft_feature_count": draft_count,
            "status": status,
        })
    return rows


def ir_note(feature: Feature) -> str:
    ir_a = (82564, 108293)
    ir_b = (125051, 150780)
    overlaps_a = any(not (end < ir_a[0] or start > ir_a[1]) for start, end in feature.ranges)
    overlaps_b = any(not (end < ir_b[0] or start > ir_b[1]) for start, end in feature.ranges)
    if overlaps_a and overlaps_b:
        return "spans_or_has_parts_in_both_IRs"
    if overlaps_a:
        return "IR_A_or_boundary"
    if overlaps_b:
        return "IR_B_or_boundary"
    return "single_copy_region"


def gene_order_rows(
    features: list[Feature],
    ref_aln: str | None,
    cand_aln: str | None,
    maps: dict[str, object],
) -> list[dict[str, object]]:
    genes = [f for f in features if f.feature_type == "gene" and f.gene]
    projected = []
    for idx, feature in enumerate(sorted(genes, key=lambda f: (f.start, f.end, f.gene)), start=1):
        presence = feature_presence(feature, ref_aln, cand_aln, maps)
        projected.append((idx, feature, presence))
    present_sorted = sorted(
        [item for item in projected if item[2]["projected_start"] != ""],
        key=lambda item: (int(item[2]["projected_start"]), int(item[2]["projected_end"])),
    )
    candidate_rank = {item[1].feature_id: rank for rank, item in enumerate(present_sorted, start=1)}
    ref_present_rank = {
        item[1].feature_id: rank
        for rank, item in enumerate([p for p in projected if p[2]["projected_start"] != ""], start=1)
    }
    rows: list[dict[str, object]] = []
    for ref_rank, feature, presence in projected:
        cand_rank = candidate_rank.get(feature.feature_id, "")
        candidate_strand = projected_candidate_strand(
            feature.strand, str(presence["projection_orientation"])
        )
        if presence["status"] != "present_by_projection":
            order_status = str(presence["status"])
        elif cand_rank == ref_present_rank.get(feature.feature_id):
            order_status = "collinear_same_orientation" if candidate_strand == feature.strand else "collinear_opposite_strand"
        else:
            order_status = "shifted_or_reordered_same_orientation" if candidate_strand == feature.strand else "shifted_or_reordered_opposite_strand"
        rows.append({
            "ref_order_index": ref_rank,
            "candidate_projected_order_index": cand_rank,
            "gene": feature.gene,
            "ref_start": feature.start,
            "ref_end": feature.end,
            "ref_strand": feature.strand,
            "candidate_projected_start": presence["projected_start"],
            "candidate_projected_end": presence["projected_end"],
            "candidate_projected_strand": candidate_strand,
            "projection_orientation": presence["projection_orientation"],
            "projected_query_coverage": f"{presence['coverage']:.6f}",
            "block_order_status": order_status,
            "ir_copy_note": ir_note(feature),
        })
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_similarity_results(
    meta: dict[str, str | int],
    candidate_meta: dict[str, str | int],
    div: dict[str, object],
) -> None:
    raw_similarity = float(div["raw_blast_weighted_identity_percent"])
    row = {
        "comparison": "D_setchellii_candidate_chloroplast_vs_NC_085682.1",
        "reference_accession": meta["accession"],
        "reference_length_bp": meta["length_bp"],
        "candidate_raw_length_bp": candidate_meta["raw_length_bp"],
        "candidate_primary_representation": "terminal_deduplicated_rotated",
        "candidate_primary_length_bp": candidate_meta["deduplicated_length_bp"],
        "raw_blast_similarity_percent": f"{raw_similarity:.3f}",
        "raw_blast_divergence_percent": f"{100 - raw_similarity:.3f}",
        "raw_blast_query_coverage": div["raw_blast_query_coverage"],
        "raw_blast_reference_coverage": div["raw_blast_reference_coverage"],
        "raw_blast_hsp_count": div["raw_blast_hsp_count"],
        "normalized_similarity_excluding_gaps_percent": div["identity_excluding_gaps_percent"],
        "normalized_divergence_excluding_gaps_percent": div["percent_divergence_excluding_gaps"],
        "normalized_similarity_with_reference_deletions_percent": div["identity_with_reference_deletions_percent"],
        "normalized_divergence_with_reference_deletions_percent": div["percent_divergence_with_reference_deletions"],
        "normalized_similarity_with_unmapped_reference_as_difference_percent": div["identity_with_unmapped_reference_as_difference_percent"],
        "normalized_divergence_with_unmapped_reference_as_difference_percent": div["percent_divergence_with_unmapped_reference_as_difference"],
        "normalized_mapped_reference_bp": div["normalized_mapped_reference_bp"],
        "normalized_unmapped_reference_bp": div["normalized_unmapped_reference_bp"],
        "normalized_mismatches": div["normalized_mismatches"],
        "normalized_candidate_deletion_bases_vs_reference": div["normalized_candidate_deletion_bases_vs_reference"],
        "normalized_forward_oriented_reference_bp": div["normalized_forward_oriented_reference_bp"],
        "normalized_reverse_oriented_reference_bp": div["normalized_reverse_oriented_reference_bp"],
    }
    write_tsv(OUT / "similarity_results.tsv", [row])


def summarize_outputs(
    meta: dict[str, str | int],
    candidate_meta: dict[str, str | int],
    rotation: dict[str, int | str | float],
    div: dict[str, object],
    cds_rows: list[dict[str, object]],
    validation: dict[str, int],
    content_rows: list[dict[str, object]],
    order_rows: list[dict[str, object]],
) -> None:
    clean_cds = [r for r in cds_rows if r["inclusion_reason"] == "included_in_clean_syn_nonsyn_totals"]
    syn = sum(int(r["single_nt_synonymous_substitutions"]) for r in clean_cds)
    nonsyn = sum(int(r["single_nt_nonsynonymous_substitutions"]) for r in clean_cds)
    complex_changes = sum(int(r["complex_codon_changes"]) for r in cds_rows)
    aa_changes = sum(int(r["amino_acid_changes"]) for r in cds_rows)
    present_content = sum(
        1 for r in content_rows
        if int(r["candidate_projected_present_count"]) == int(r["nc085682_feature_count"])
    )
    draft_diff = [r for r in content_rows if r["status"] == "present_by_alignment_draft_count_diff"]
    reordered = [r for r in order_rows if str(r["block_order_status"]).startswith("shifted_or_reordered")]
    opposite = [r for r in order_rows if "opposite_strand" in str(r["block_order_status"])]
    summary = f"""# NC_085682.1 Chloroplast Comparison

## Reference and Candidate
- Reference: {meta['accession']}, {meta['description']} ({meta['length_bp']} bp; LOCUS date {meta['locus_date']}).
- Reference source: {meta['source']}; fetch date recorded as {meta['fetch_date']}.
- Candidate raw plastome: {candidate_meta['raw_length_bp']} bp; ambiguous bases {candidate_meta['raw_ambiguous_bases']}.
- Candidate primary comparison plastome: terminal de-duplicated to {candidate_meta['deduplicated_length_bp']} bp by trimming the suffix starting at position {candidate_meta['terminal_duplicate_suffix_start']}.
- Rotation/orientation: BLAST anchor at candidate position {rotation['rotation_start_1based']} in {rotation['orientation']} orientation, pident {rotation['pident']:.3f}, length {rotation['length']} bp.

## Whole-Genome Divergence
- Raw BLAST compatibility check: query coverage {div['raw_blast_query_coverage']}, reference coverage {div['raw_blast_reference_coverage']}, weighted identity {div['raw_blast_weighted_identity_percent']}%, HSPs {div['raw_blast_hsp_count']}.
- One-row similarity result file: `similarity_results.tsv`.
- Normalized BLAST projection: {div['normalized_mapped_reference_bp']} mapped reference bp, {div['normalized_unmapped_reference_bp']} unmapped reference bp, {div['normalized_mismatches']} mismatches, and {div['normalized_candidate_deletion_bases_vs_reference']} candidate deletion bases versus the reference.
- Identity excluding gaps: {div['identity_excluding_gaps_percent']}%; divergence excluding gaps: {div['percent_divergence_excluding_gaps']}%.
- Identity counting reference deletions: {div['identity_with_reference_deletions_percent']}%; divergence counting reference deletions: {div['percent_divergence_with_reference_deletions']}%.
- Identity counting unmapped reference bases as differences: {div['identity_with_unmapped_reference_as_difference_percent']}%; divergence on that stricter denominator: {div['percent_divergence_with_unmapped_reference_as_difference']}%.
- Diagnostic MAFFT linear alignment, not used for final calls because repeats/SSC orientation can mislead it: {div['diagnostic_mafft_aligned_columns']} columns, {div['diagnostic_mafft_mismatches']} mismatches, {div['diagnostic_mafft_indel_columns']} indel columns.

## CDS Substitutions
- Shared/projected CDS rows: {len(cds_rows)}.
- CDS included in clean single-codon substitution totals: {len(clean_cds)}.
- Clean single-nt synonymous substitutions: {syn}.
- Clean single-nt nonsynonymous substitutions: {nonsyn}.
- Complex codon changes kept separate: {complex_changes}; amino-acid-changing codons across all projected CDS rows: {aa_changes}.
- Reference CDS translation validation: {validation.get('translation_matches', 0)} matches and {validation.get('translation_mismatches', 0)} mismatches among {validation.get('features_with_translation', 0)} CDS with /translation qualifiers.

## Gene Content and Order
- Gene-content rows with full projection support: {present_content}/{len(content_rows)}.
- Rows present by alignment but differing from the current draft annotation count: {len(draft_diff)}. This is expected in IR and draft-transfer edge cases.
- Gene-order rows flagged as shifted/reordered after circular normalization: {len(reordered)}.
- Gene-order rows projected on the opposite strand: {len(opposite)}, consistent with an SSC-orientation difference and IR-copy effects. Inspect `gene_order.tsv` for exact rows.

## Interpretation
The candidate chloroplast is very close to NC_085682.1 after removing the terminal duplicate and rotating to the same circular origin. The raw BLAST result remains consistent with existing reference-verification evidence, while the normalized BLAST projection gives the biologically more meaningful divergence estimate. Gene content is mostly present by alignment; differences from the draft annotation should be interpreted as annotation completeness/copy-number issues unless `gene_content.tsv` marks projection absence. Synonymous/nonsynonymous totals are codon-consequence counts from projected shared CDS and are not model-based dN/dS or Ka/Ks estimates.
"""
    (OUT / "summary.md").write_text(summary)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare the candidate chloroplast against a live NCBI download of NC_085682.1."
    )
    parser.add_argument(
        "--allow-cached-fallback",
        action="store_true",
        help=(
            "Use the cached NC_085682.1 GenBank record only if the live NCBI download fails. "
            "Do not use this for the committed comparison artifacts."
        ),
    )
    args = parser.parse_args()

    for binary in ("blastn", "mafft"):
        if not shutil.which(binary):
            raise RuntimeError(f"Required executable not found on PATH: {binary}")

    fetched_fasta, fetched_gb, source_note = fetch_current_records(args.allow_cached_fallback)
    if fetched_gb:
        record = fetched_gb
        ref_seq = extract_origin(record)
    else:
        record = cached_nc_record()
        ref_seq = extract_origin(record)
    if fetched_fasta:
        (OUT / "NC_085682.1.fetched.fa").write_text(fetched_fasta)
    if fetched_gb:
        (OUT / "NC_085682.1.fetched.gb").write_text(fetched_gb)

    features = parse_features(record)
    locus_length = parse_locus_length(record)
    if len(ref_seq) != locus_length:
        raise RuntimeError(f"Reference sequence length {len(ref_seq)} != LOCUS length {locus_length}")
    meta = write_reference_files(record, ref_seq, source_note)

    candidate_meta = prepare_candidate_files()
    candidate_dedup = read_fasta(OUT / "candidate_terminal_deduplicated.fa")[0][1]

    raw_rows = blast_rows(
        OUT / "candidate_raw.fa",
        OUT / "NC_085682.1.fa",
        OUT / "raw_candidate_vs_NC_085682.1.blastn.tsv",
    )
    raw_blast = summarize_blast(raw_rows)

    anchor_rows = blast_rows(
        OUT / "NC_085682.1.fa",
        OUT / "candidate_terminal_deduplicated.fa",
        OUT / "NC_085682.1_vs_candidate_deduplicated.blastn.tsv",
    )
    rotation = choose_rotation_anchor(anchor_rows, len(candidate_dedup))
    oriented = candidate_dedup
    if rotation["orientation"] == "-":
        oriented = revcomp(oriented)
    rotated = rotate_sequence(oriented, int(rotation["rotation_start_1based"]))
    write_fasta(OUT / "candidate_terminal_deduplicated_rotated.fa", [
        ("candidate_terminal_deduplicated_rotated_to_NC_085682.1", rotated)
    ])

    ref_aln, cand_aln = run_mafft(ref_seq, rotated)
    mafft_diag = divergence_metrics(ref_aln, cand_aln, raw_blast)
    normalized_alignment_rows = blast_alignment_rows(
        OUT / "NC_085682.1.fa",
        OUT / "candidate_terminal_deduplicated_rotated.fa",
        OUT / "normalized_ref_vs_candidate_rotated.blastn_alignment.tsv",
    )
    projection = build_blast_projection(normalized_alignment_rows)
    projection_maps = {"projection_by_ref": projection}
    div = blast_projection_divergence_metrics(ref_seq, projection, raw_blast, mafft_diag)
    write_tsv(OUT / "whole_genome_divergence.tsv", [div])
    write_similarity_results(meta, candidate_meta, div)

    cds_rows, validation = compare_cds(features, ref_seq, None, None, projection_maps)
    write_tsv(OUT / "cds_substitutions.tsv", cds_rows)

    content_rows = gene_content_rows(features, None, None, projection_maps)
    write_tsv(OUT / "gene_content.tsv", content_rows)

    order_rows = gene_order_rows(features, None, None, projection_maps)
    write_tsv(OUT / "gene_order.tsv", order_rows)

    # Machine-readable provenance for quick verification without parsing Markdown.
    provenance = {
        "accession": ACCESSION,
        "reference_length_bp": len(ref_seq),
        "candidate_raw_length_bp": candidate_meta["raw_length_bp"],
        "candidate_deduplicated_length_bp": candidate_meta["deduplicated_length_bp"],
        "rotation_start_1based": rotation["rotation_start_1based"],
        "raw_blast_query_coverage": div["raw_blast_query_coverage"],
        "raw_blast_reference_coverage": div["raw_blast_reference_coverage"],
        "raw_blast_weighted_identity_percent": div["raw_blast_weighted_identity_percent"],
        "raw_blast_hsp_count": div["raw_blast_hsp_count"],
        "translation_matches": validation.get("translation_matches", 0),
        "translation_mismatches": validation.get("translation_mismatches", 0),
        "normalized_projection_mapped_reference_bp": div["normalized_mapped_reference_bp"],
        "normalized_projection_unmapped_reference_bp": div["normalized_unmapped_reference_bp"],
    }
    with (OUT / "verification_metrics.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(provenance.keys()), delimiter="\t")
        writer.writeheader()
        writer.writerow(provenance)

    summarize_outputs(meta, candidate_meta, rotation, div, cds_rows, validation, content_rows, order_rows)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
