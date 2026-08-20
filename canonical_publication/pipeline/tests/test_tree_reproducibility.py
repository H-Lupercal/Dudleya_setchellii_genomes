from pathlib import Path

from organelle_pipeline.tree_reproducibility import compare_unrooted_trees


def write_tree(path: Path, newick: str) -> Path:
    path.write_text(newick + "\n")
    return path


def test_equivalent_unrooted_orderings_reproduce(tmp_path: Path) -> None:
    left = write_tree(tmp_path / "left.tree", "((A:1,B:1)90/100:1,(C:1,D:1)90/100:1);")
    right = write_tree(tmp_path / "right.tree", "((D:1,C:1)90/100:1,(B:1,A:1)90/100:1);")
    comparison = compare_unrooted_trees(left, right)
    assert comparison.full_unrooted_rf == 0
    assert comparison.strong_split_symmetric_difference == 0
    assert comparison.strong_topology_reproduced


def test_weak_polytomy_resolutions_are_reported_but_do_not_fail(tmp_path: Path) -> None:
    left = write_tree(tmp_path / "left.tree", "(((A:1,B:1)0/0:0,C:1),D:1,E:1);")
    right = write_tree(tmp_path / "right.tree", "(((A:1,C:1)0/0:0,B:1),D:1,E:1);")
    comparison = compare_unrooted_trees(left, right)
    assert comparison.full_unrooted_rf > 0
    assert comparison.strong_split_symmetric_difference == 0
    assert comparison.strong_topology_reproduced


def test_different_strong_splits_fail_reproducibility(tmp_path: Path) -> None:
    left = write_tree(tmp_path / "left.tree", "(((A:1,B:1)90/100:1,C:1),D:1,E:1);")
    right = write_tree(tmp_path / "right.tree", "(((A:1,C:1)90/100:1,B:1),D:1,E:1);")
    comparison = compare_unrooted_trees(left, right)
    assert comparison.strong_split_symmetric_difference > 0
    assert not comparison.strong_topology_reproduced
