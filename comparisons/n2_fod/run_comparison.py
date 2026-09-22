"""Generate an N2 N_FOD comparison between ORCA and PySCF.

ORCA is run with its ``! FOD`` defaults.  PySCF uses the corresponding
TPSS/def2-TZVP calculation at 5000 K.  Program-specific numerical defaults
are intentionally left unchanged.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import re
import shutil
import subprocess

from pyscf import gto

from static_correlation_thermo import run_ftdft


DISTANCES_ANGSTROM = (1.10, 1.40, 1.80, 2.20, 2.60, 3.00)
N_FOD_PATTERN = re.compile(
    r"\bN_FOD\s*=\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)"
)


def orca_input(distance: float) -> str:
    """Return a minimal ORCA input using the program's FOD defaults."""
    half_distance = distance / 2.0
    return (
        "! FOD\n\n"
        "* xyz 0 1\n"
        f"N  0.0  0.0  {-half_distance:.8f}\n"
        f"N  0.0  0.0   {half_distance:.8f}\n"
        "*\n"
    )


def parse_orca_n_fod(output_file: Path) -> float:
    """Extract the final N_FOD value from an ORCA output file."""
    matches = N_FOD_PATTERN.findall(output_file.read_text(encoding="utf-8"))
    if not matches:
        raise RuntimeError(f"N_FOD was not found in {output_file}")
    return float(matches[-1])


def run_orca(orca: Path, distance: float, run_directory: Path) -> float:
    """Run one default ORCA FOD calculation and return N_FOD."""
    run_directory.mkdir(parents=True, exist_ok=True)
    input_file = run_directory / "n2.inp"
    output_file = run_directory / "n2.out"
    input_file.write_text(orca_input(distance), encoding="utf-8")
    with output_file.open("w", encoding="utf-8") as output:
        completed = subprocess.run(
            [str(orca.resolve()), input_file.name],
            cwd=run_directory,
            stdout=output,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if completed.returncode != 0:
        raise RuntimeError(
            f"ORCA failed for R={distance:.2f} A; inspect {output_file}"
        )
    return parse_orca_n_fod(output_file)


def run_pyscf(distance: float) -> float:
    """Run the PySCF calculation matching ORCA's default FOD model chemistry."""
    half_distance = distance / 2.0
    molecule = gto.M(
        atom=(
            f"N 0 0 {-half_distance:.8f}; "
            f"N 0 0 {half_distance:.8f}"
        ),
        basis="def2-tzvp",
        unit="Angstrom",
        spin=0,
        verbose=0,
    )
    result = run_ftdft(
        molecule,
        xc="tpss",
        temperature=5000.0,
    )
    return result.n_fod


def find_orca(explicit_path: Path | None) -> Path:
    if explicit_path is not None:
        return explicit_path
    executable = shutil.which("orca")
    if executable is None:
        raise RuntimeError("ORCA is not on PATH; supply --orca /path/to/orca")
    return Path(executable)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orca", type=Path, help="path to the ORCA executable")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("n2_fod_runs"),
        help="directory for ORCA runs and the comparison CSV",
    )
    args = parser.parse_args()
    orca = find_orca(args.orca)

    rows: list[tuple[float, float, float, float, float]] = []
    for distance in DISTANCES_ANGSTROM:
        point_directory = args.output_dir / f"n2_{distance:.2f}A"
        orca_n_fod = run_orca(orca, distance, point_directory)
        pyscf_n_fod = run_pyscf(distance)
        difference = pyscf_n_fod - orca_n_fod
        relative_difference = difference / orca_n_fod if orca_n_fod else float("nan")
        rows.append(
            (distance, orca_n_fod, pyscf_n_fod, difference, relative_difference)
        )
        print(
            f"R={distance:.2f} A  ORCA={orca_n_fod:.6f}  "
            f"PySCF={pyscf_n_fod:.6f}  delta={difference:+.6f}"
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_file = args.output_dir / "n2_fod_comparison.csv"
    with csv_file.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "distance_angstrom",
                "orca_n_fod",
                "pyscf_n_fod",
                "pyscf_minus_orca",
                "relative_difference",
            ]
        )
        writer.writerows(rows)
    print(f"Wrote {csv_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

