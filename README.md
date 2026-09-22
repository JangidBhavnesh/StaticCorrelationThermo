# Finite-temperature DFT in PySCF

This small package runs molecular Kohn--Sham DFT with Fermi--Dirac orbital
occupations.  It is intended to make near-degeneracy (static-correlation)
effects visible for a chosen geometry and exchange-correlation functional.

At an electronic temperature `T`, PySCF minimizes an approximate electronic
Helmholtz free energy

```text
A(T) = E(T) - k_B T S
```

and reports fractional KS occupations and the noninteracting entropy `S`.
When a bond is stretched and its bonding and antibonding orbitals become nearly
degenerate, both orbitals acquire appreciable fractional occupation.  This is a
useful qualitative static-correlation diagnostic.  The temperature is a
fictitious regularization parameter here, not a nuclear temperature.

## Install and run

```bash
python -m pip install -e .
ftdft examples/h2_stretched.xyz --xc pbe --basis def2-svp \
  --temperature 0 5000
```

The table contains the KS internal energy `E(T)`, electronic free energy
`A(T)`, entropy in units of `k_B`, and a fractional-occupation index.  The
index is zero for integer restricted occupations and approaches 2 when two
spatial orbitals each have occupation 1.  Fractionally occupied orbitals are
listed below each finite-temperature result.

To see the effect develop as H2 is stretched, run:

```bash
python examples/h2_dissociation.py
```

The Python interface accepts any built PySCF molecule and any XC expression
understood by PySCF:

```python
from pyscf import gto
from static_correlation_thermo import run_ftdft

mol = gto.M(atom="H 0 0 -1.5; H 0 0 1.5", basis="def2-svp")
result = run_ftdft(mol, xc="pbe0", temperature=5000.0)

print(result.free_energy, result.entropy_kb)
print(result.fractional_orbitals())
```

For open-shell systems the helper selects UKS and conserves the alpha and beta
electron counts separately.  Use a spin-restricted calculation for stretched
closed-shell H2 if the goal is to expose static correlation; unrestricted DFT
can instead lower its energy by breaking spin symmetry.

## Important limitation

This is a Mermin-style finite-temperature KS calculation using PySCF's
smearing support, but standard PySCF functionals are ground-state XC
approximations.  Consequently it is **not** a fully temperature-dependent XC
theory, and fractional KS occupations are not correlated natural-orbital
occupations.  Treat the entropy and fractional-occupation index as qualitative
diagnostics and check their dependence on the chosen fictitious temperature,
functional, basis, and spin constraint.

Run the tests with `pytest`.
