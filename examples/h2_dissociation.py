"""Show the growth and location of static correlation as H2 dissociates."""

from pathlib import Path

from pyscf import gto

from static_correlation_thermo import run_ftdft


output_directory = Path("fod_results")

print(" R / A       A(T) / Eh       S / kB       N_FOD")
for distance in (0.74, 1.50, 2.00, 3.00, 4.00):
    molecule = gto.M(
        atom=f"H 0 0 0; H 0 0 {distance}",
        basis="def2-svp",
        unit="Angstrom",
        verbose=0,
    )
    result = run_ftdft(
        molecule,
        xc="pbe",
        temperature=5000.0,
        output_prefix=output_directory / f"h2_{distance:.2f}A",
    )
    print(
        f"{distance:6.2f}  {result.free_energy:15.8f}  "
        f"{result.entropy_kb:11.6f}  {result.n_fod:11.6f}"
    )

print(f"\nMolden and FOD cube files were written to {output_directory}/")
