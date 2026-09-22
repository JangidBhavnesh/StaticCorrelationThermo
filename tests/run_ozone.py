from pyscf import gto
from static_correlation_thermo import run_ftdft

mol = gto.M(
    atom="""
    O   0.00000000017911   0.00000000000000   0.43702029765776
    O  -1.09512651993192   0.00000000000000  -0.21851064888325
    O   1.09512651975281   0.00000000000000  -0.21851064877451
    """,
    unit="Angstrom",
    charge=0,
    spin=0,
    basis="def2-tzvp"
)

result = run_ftdft(
    mol,
    xc="tpss",
    temperature=5000.0,
    output_prefix="ozone",
)

print(result.free_energy)
print(result.entropy_kb)
print(result.n_fod)
print(result.molden_file)
print(result.fod_cube_file)
