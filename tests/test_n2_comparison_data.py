import csv
from pathlib import Path

import pytest


def test_committed_n2_comparison_data_are_self_consistent():
    data_file = (
        Path(__file__).resolve().parents[1]
        / "comparisons"
        / "n2_fod"
        / "n2_fod_comparison.csv"
    )
    with data_file.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 6
    differences = []
    for row in rows:
        orca = float(row["orca_n_fod"])
        pyscf = float(row["pyscf_n_fod"])
        difference = float(row["pyscf_minus_orca"])
        assert difference == pytest.approx(pyscf - orca, abs=1.0e-15)
        differences.append(abs(difference))

    assert max(differences) < 1.0e-3

