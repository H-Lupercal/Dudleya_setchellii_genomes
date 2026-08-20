from organelle_pipeline.annotations import gff_phase, projected_interval
from organelle_pipeline.reference_evidence import BlastHit


def test_projected_interval_preserves_reverse_orientation() -> None:
    hit = BlastHit("feature", "cp", 99.0, 100, 1, 100, 500, 401, 200, 100, 1000)
    assert projected_interval(hit) == (400, 500, "-")


def test_cds_codon_start_is_converted_to_gff3_phase() -> None:
    assert gff_phase("CDS", "1") == "0"
    assert gff_phase("CDS", "2") == "1"
    assert gff_phase("CDS", "3") == "2"
    assert gff_phase("gene", None) == "."
