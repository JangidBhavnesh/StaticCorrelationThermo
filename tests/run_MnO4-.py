from pyscf import gto
from static_correlation_thermo import run_ftdft

mol = gto.M(
    atom="""
    Mn   0.00000   0.00000   0.00000
    O    0.93250  -0.93250   0.93250
    O   -0.93250   0.93250   0.93250
    O   -0.93250  -0.93250  -0.93250
    O    0.93250   0.93250  -0.93250
    """,
    unit="Angstrom",
    charge=-1,
    spin=0,
    basis="def2-tzvp"
)

result = run_ftdft(
    mol,
    xc="tpss",
    temperature=5000.0,
    output_prefix="permanganate",
)

print(result.free_energy)
print(result.entropy_kb)
print(result.n_fod)
print(result.molden_file)
print(result.fod_cube_file)
