from organelle_pipeline.qc import (
    SampleBreadth,
    build_depth_command,
    select_organelle_samples,
    summarize_depths,
    summarize_fastp_report,
    summarize_masked_depths,
    support_intervals,
)


def test_fastp_summary_validates_and_exposes_read_filtering_metrics() -> None:
    observed = summarize_fastp_report(
        {
            "summary": {
                "before_filtering": {
                    "total_reads": 100,
                    "total_bases": 10_000,
                    "q20_rate": 0.90,
                    "q30_rate": 0.80,
                },
                "after_filtering": {
                    "total_reads": 80,
                    "total_bases": 7_500,
                    "q20_rate": 0.98,
                    "q30_rate": 0.92,
                },
            },
            "filtering_result": {
                "passed_filter_reads": 80,
                "low_quality_reads": 10,
                "too_many_N_reads": 2,
                "adapter_dimer_reads": 3,
                "too_short_reads": 5,
            },
            "adapter_cutting": {
                "adapter_trimmed_reads": 25,
                "adapter_trimmed_bases": 200,
            },
            "duplication": {"rate": 0.1},
        }
    )

    assert observed.input_reads == 100
    assert observed.passing_reads == 80
    assert observed.read_retention == 0.8
    assert observed.adapter_trimmed_reads == 25
    assert observed.duplication_rate == 0.1


def test_organelle_eligibility_is_independent_and_shared_is_intersection() -> None:
    rows = [
        SampleBreadth("both", cp_dp5=0.91, mt_dp5=0.88),
        SampleBreadth("cp_only", cp_dp5=0.80, mt_dp5=0.79),
        SampleBreadth("mt_only", cp_dp5=0.50, mt_dp5=0.95),
    ]

    selected = select_organelle_samples(rows, minimum_dp5_breadth=0.80)

    assert selected.cp == ("both", "cp_only")
    assert selected.mt == ("both", "mt_only")
    assert selected.shared == ("both",)


def test_threshold_is_inclusive() -> None:
    selected = select_organelle_samples(
        [SampleBreadth("edge", cp_dp5=0.80, mt_dp5=0.80)],
        minimum_dp5_breadth=0.80,
    )
    assert selected.shared == ("edge",)


def test_depth_summary_reports_all_requested_breadths() -> None:
    summary = summarize_depths([0, 1, 3, 5, 10], thresholds=(1, 3, 5, 10))
    assert summary.mean_depth == 3.8
    assert summary.breadth == {1: 0.8, 3: 0.6, 5: 0.4, 10: 0.2}


def test_analysis_mask_breadth_excludes_structurally_unassignable_repeat_sites() -> None:
    summary = summarize_masked_depths([5, 0, 5], [True, False, True], thresholds=(5,))
    assert summary.reference_length == 2
    assert summary.breadth == {5: 1.0}


def test_samtools_depth_quality_flags_are_not_reversed() -> None:
    command = build_depth_command(
        "sample.bam",
        "chloroplast",
        minimum_mapping_quality=27,
        minimum_base_quality=13,
    )

    assert command[command.index("-q") + 1] == "13"
    assert command[command.index("-Q") + 1] == "27"
    assert "-s" in command


def test_support_intervals_use_inclusive_sample_fraction_and_minimum_length() -> None:
    intervals = support_intervals(
        [8, 8, 7, 8, 8, 8, 0, 8, 8],
        sample_count=10,
        minimum_fraction=0.8,
        minimum_length=2,
    )
    assert intervals == ((0, 2), (3, 6), (7, 9))
