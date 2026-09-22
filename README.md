# Static correlation Thermometer with FT-DFT in PySCF

This project will provide a small PySCF workflow for visualizing static
correlation at a supplied molecular geometry and exchange-correlation (XC)
functional. The calculation uses Fermi--Dirac occupations at a **fictitious**
electronic temperature, making near-degenerate frontier orbitals fractionally
occupied.

The repository currently contains the package scaffold and design described
below. The FT-DFT calculation, FOD-density export, and examples are the next
implementation steps.

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

For orbital energy $\epsilon_i$, chemical potential $\mu$, and
$\sigma = k_B T$, the occupation of one spin orbital is


$$f_i = \frac{1}{\exp((\epsilon_i - \mu) / \sigma) + 1}$$


The electron number determines $\mu$. PySCF then reports the self-consistent
internal energy $E(T)$, dimensionless noninteracting entropy $S/k_B$, and
electronic Helmholtz free energy


$$A(T) = E(T) - \sigma S/k_B$$

As a bond is stretched, its bonding and antibonding orbitals become nearly
degenerate. In a spin-restricted calculation they approach occupations of one
electron each, revealing the onset of static correlation.

### Fractional-occupation density

Following the FOD analysis used in ORCA, the real-space fractional-occupation
density will be constructed from the finite-temperature molecular spin
orbitals:


$$\rho_\text{FOD}(r) = \sum_i w_i |\phi_i(r)|^2 \\
w_i = \begin{cases}
    1 - f_i &   \text{for } \epsilon_i < \mu,  \\
    f_i     &   \text{for } \epsilon_i > \mu
\end{cases}
    $$


Here, `f_i` is a spin-orbital occupation between zero and one. Thus, fully
occupied and fully empty orbitals make no contribution, while fractionally
occupied frontier orbitals reveal where the statically correlated ("hot")
electrons are localized. Its spatial integral gives the scalar diagnostic

$$
N_\text{FOD} = \int \rho_\text{FOD}(r) dr = \sum_i w_i
$$

This follows the definition in the
[ORCA FOD documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/fod.html)
and the original
[Grimme--Hansen FOD work](https://doi.org/10.1002/anie.201501887).

The primary outputs will therefore be:

- orbital energies and fractional occupations;
- $E(T)$, $A(T)$, and $S/k_B$;
- the integrated static-correlation diagnostic $N_\text{FOD}$;
- a Molden file containing the converged MOs, energies, spins, and fractional
  occupations;
- a Gaussian cube file containing $\rho_\text{FOD}(r)$ for direct visualization;
- an optional scan over temperature or molecular geometry.

## Visualization files

Each calculation will write two complementary files using the chosen output
prefix:

```text
<prefix>.molden     geometry, basis, MOs, energies, spins, and occupations
<prefix>.fod.cube   three-dimensional fractional-occupation density
```

The Molden file can be opened in a Molden-compatible viewer to inspect the
individual fractionally occupied frontier orbitals. The FOD itself is a sum of
weighted orbital densities, not a single molecular orbital, so it cannot be
represented faithfully as one Molden orbital. The ORCA-like FOD picture is
therefore obtained by opening `<prefix>.fod.cube` in a cube-compatible viewer
such as Molden, Jmol, Avogadro, VMD, or Chemcraft and displaying a positive
isosurface. An isovalue of `0.005 e/bohr^3` is a useful ORCA-compatible starting
point; the value may need adjustment for a particular molecule.

Internally, the implementation will form the FOD one-particle density matrix
from the weighted MO coefficients and use PySCF to evaluate it on the cube
grid. The numerical integral of the cube will be checked against $N_\text{FOD}$.

## Planned Python interface

```python
from pyscf import gto
from static_correlation_thermo import run_ftdft

mol = gto.M(
    atom="H 0 0 -1.5; H 0 0 1.5",
    basis="def2-svp",
    unit="Angstrom",
)

result = run_ftdft(
    mol,
    xc="pbe",
    temperature=5000.0,
    output_prefix="h2_stretched",
)
print(result.free_energy)
print(result.entropy_kb)
print(result.n_fod)
print(result.molden_file)   # h2_stretched.molden
print(result.fod_cube_file) # h2_stretched.fod.cube
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
- `N_FOD` and the FOD isosurface are qualitative static-correlation
  indicators, not a replacement for a multireference calculation;
- results should be checked against temperature, XC functional, basis set,
  geometry, and the restricted/unrestricted choice.

The first implementation will target molecular RKS/UKS calculations and use
PySCF's `scf.addons.smearing(..., method="fermi")` machinery.
