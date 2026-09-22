# Static correlation Thermometer with FT-DFT in PySCF

This project will provide a small PySCF workflow for visualizing static
correlation at a supplied molecular geometry and exchange-correlation (XC)
functional. The calculation uses Fermi--Dirac occupations at a **fictitious**
electronic temperature, making near-degenerate frontier orbitals fractionally
occupied.

The repository currently contains the package scaffold and design described
below. The FT-DFT calculation code and examples are the next implementation
steps.

## Repository layout

```text
StaticCorrelationThermo/
├── src/
│   └── static_correlation_thermo/
│       └── __init__.py
├── examples/
├── tests/
└── README.md
```

The importable Python package is `static_correlation_thermo`. Calculation code
will live under this package, example calculations under `examples/`, and
automated checks under `tests/`.

## Method

For orbital energy `epsilon_i`, chemical potential `mu`, and
`sigma = k_B T`, the occupation of one spin orbital is

```text
f_i = 1 / (exp((epsilon_i - mu) / sigma) + 1).
```

The electron number determines `mu`. PySCF then reports the self-consistent
internal energy `E(T)`, dimensionless noninteracting entropy `S/k_B`, and
electronic Helmholtz free energy

```text
A(T) = E(T) - sigma (S/k_B).
```

As a bond is stretched, its bonding and antibonding orbitals become nearly
degenerate. In a spin-restricted calculation they approach occupations of one
electron each, revealing the onset of static correlation. The primary outputs
will therefore be:

- orbital energies and fractional occupations;
- `E(T)`, `A(T)`, and `S/k_B`;
- a simple fractional-occupation diagnostic;
- an optional scan over temperature or molecular geometry.

## Planned Python interface

```python
from pyscf import gto
from static_correlation_thermo import run_ftdft

mol = gto.M(
    atom="H 0 0 -1.5; H 0 0 1.5",
    basis="def2-svp",
    unit="Angstrom",
)

result = run_ftdft(mol, xc="pbe", temperature=5000.0)
print(result.free_energy)
print(result.entropy_kb)
print(result.fractional_orbitals())
```

The initial example will scan the H--H distance. Near equilibrium the bonding
orbital is close to doubly occupied and the antibonding orbital is nearly
empty. At stretched geometries both occupations move toward one as the
bonding--antibonding gap closes.

## Spin choice

Closed-shell examples intended to expose static correlation should remain
spin restricted. Unrestricted DFT can instead lower the energy by breaking
spin symmetry and may hide the fractional-occupation signature. Open-shell
systems will use UKS with the alpha and beta electron counts fixed separately.

## Scope and limitations

This workflow is a Mermin-style finite-temperature KS calculation based on
PySCF's Fermi smearing support. Standard PySCF XC functionals are ground-state
approximations; they do not acquire explicit temperature dependence here.
Accordingly:

- the electronic temperature is a diagnostic/regularization parameter, not
  the nuclear temperature of an experiment;
- fractional KS occupations are not correlated natural-orbital occupations;
- entropy and fractional occupations are qualitative static-correlation
  indicators, not a replacement for a multireference calculation;
- results should be checked against temperature, XC functional, basis set,
  geometry, and the restricted/unrestricted choice.

The first implementation will target molecular RKS/UKS calculations and use
PySCF's `scf.addons.smearing(..., method="fermi")` machinery.
