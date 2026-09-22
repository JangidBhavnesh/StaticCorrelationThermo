"""Finite-temperature molecular DFT and FOD analysis using PySCF.

PySCF supplies the Kohn--Sham solver, Fermi--Dirac smearing, Molden writer,
and cube-grid density evaluator.  This module adds the fractional-occupation
density (FOD) weights of Grimme and Hansen and a small result-oriented API.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from typing import Any
import warnings

import numpy as np
from numpy.typing import ArrayLike, NDArray
from pyscf import dft, gto, scf
from pyscf.data import nist
from pyscf.tools import cubegen, molden


BOLTZMANN_HARTREE_PER_KELVIN = nist.BOLTZMANN / nist.HARTREE2J
"""Boltzmann constant in hartree per kelvin, using PySCF's constants."""


def kelvin_to_hartree(temperature: float) -> float:
    """Convert an electronic temperature in kelvin to ``k_B T`` in hartree."""
    temperature = float(temperature)
    if not math.isfinite(temperature) or temperature < 0.0:
        raise ValueError("temperature must be a finite, non-negative number")
    return temperature * BOLTZMANN_HARTREE_PER_KELVIN


def fod_weights(mo_occ: ArrayLike, *, unrestricted: bool = False) -> NDArray[np.float64]:
    """Return ORCA-style FOD weights from PySCF orbital occupations.

    Restricted PySCF occupations lie in [0, 2] and represent two equivalent
    spin orbitals.  Their combined FOD weight is ``min(n, 2-n)``.
    Unrestricted occupations lie in [0, 1] for each spin and have weight
    ``min(f, 1-f)``.  These expressions are equivalent to using holes below
    the chemical potential and particles above it for Fermi occupations.
    """
    occupations = np.asarray(mo_occ, dtype=float)
    upper = 1.0 if unrestricted else 2.0
    expected_ndim = 2 if unrestricted else 1
    if occupations.ndim != expected_ndim:
        kind = "two-dimensional" if unrestricted else "one-dimensional"
        raise ValueError(f"occupations must be {kind}")
    tolerance = 1.0e-10
    if np.any(occupations < -tolerance) or np.any(occupations > upper + tolerance):
        raise ValueError(f"occupations must lie between 0 and {upper:g}")
    occupations = np.clip(occupations, 0.0, upper)
    return np.minimum(occupations, upper - occupations)


def fod_density_matrix(
    mo_coeff: ArrayLike,
    mo_occ: ArrayLike,
    *,
    unrestricted: bool = False,
) -> NDArray[Any]:
    """Construct the AO-basis FOD density matrix from PySCF MOs.

    For UKS, the alpha and beta contributions are added because the exported
    cube contains the total FOD density.
    """
    coefficients = np.asarray(mo_coeff)
    weights = fod_weights(mo_occ, unrestricted=unrestricted)

    if unrestricted:
        if coefficients.ndim != 3 or coefficients.shape[0] != 2:
            raise ValueError("unrestricted MO coefficients must have shape (2, nao, nmo)")
        if coefficients.shape[2] != weights.shape[1]:
            raise ValueError("MO coefficients and occupations are inconsistent")
        density_matrix = sum(
            (coefficients[spin] * weights[spin]) @ coefficients[spin].conj().T
            for spin in range(2)
        )
    else:
        if coefficients.ndim != 2:
            raise ValueError("restricted MO coefficients must have shape (nao, nmo)")
        if coefficients.shape[1] != weights.shape[0]:
            raise ValueError("MO coefficients and occupations are inconsistent")
        density_matrix = (coefficients * weights) @ coefficients.conj().T

    density_matrix = np.real_if_close(density_matrix)
    if np.iscomplexobj(density_matrix):
        raise ValueError("complex-valued FOD densities are not supported by cube export")
    return np.asarray(density_matrix, dtype=float)


@dataclass(frozen=True)
class FTDFTResult:
    """Observable values and output paths from one FT-DFT calculation."""

    temperature_kelvin: float
    sigma_hartree: float
    internal_energy: float
    free_energy: float
    zero_temperature_energy: float
    entropy_kb: float
    n_fod: float
    cube_integral: float | None
    mo_energy: NDArray[np.float64]
    mo_occ: NDArray[np.float64]
    unrestricted: bool
    molden_file: Path | None
    fod_cube_file: Path | None
    mean_field: Any

    def fractional_orbitals(
        self, cutoff: float = 1.0e-4
    ) -> list[tuple[str, int, float, float, float]]:
        """List fractionally occupied orbitals.

        Each row is ``(spin, zero_based_index, energy, occupation, FOD_weight)``.
        """
        if not 0.0 < cutoff < 0.5:
            raise ValueError("cutoff must lie between 0 and 0.5")
        weights = fod_weights(self.mo_occ, unrestricted=self.unrestricted)
        rows: list[tuple[str, int, float, float, float]] = []
        if self.unrestricted:
            for spin_index, spin_name in enumerate(("alpha", "beta")):
                for index, (energy, occupation, weight) in enumerate(
                    zip(
                        self.mo_energy[spin_index],
                        self.mo_occ[spin_index],
                        weights[spin_index],
                        strict=True,
                    )
                ):
                    if weight > cutoff:
                        rows.append(
                            (spin_name, index, float(energy), float(occupation), float(weight))
                        )
        else:
            for index, (energy, occupation, weight) in enumerate(
                zip(self.mo_energy, self.mo_occ, weights, strict=True)
            ):
                if weight > cutoff:
                    rows.append(
                        ("restricted", index, float(energy), float(occupation), float(weight))
                    )
        return rows


def _cube_integral(
    mol: gto.Mole,
    density: NDArray[np.float64],
    cube_points: tuple[int, int, int],
    resolution: float | None,
    margin: float,
) -> float:
    """Integrate a PySCF cube grid with the trapezoidal rule."""
    cube = cubegen.Cube(
        mol,
        nx=cube_points[0],
        ny=cube_points[1],
        nz=cube_points[2],
        resolution=resolution,
        margin=margin,
    )
    integral = np.trapz(density, x=cube.zs, axis=2)
    integral = np.trapz(integral, x=cube.ys, axis=1)
    integral = np.trapz(integral, x=cube.xs, axis=0)
    return float(integral * abs(np.linalg.det(cube.box)))


def _output_paths(output_prefix: str | Path) -> tuple[Path, Path]:
    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    return Path(f"{prefix}.molden"), Path(f"{prefix}.fod.cube")


def run_ftdft(
    mol: gto.Mole,
    *,
    xc: str = "pbe",
    temperature: float = 5000.0,
    output_prefix: str | Path | None = None,
    unrestricted: bool | None = None,
    grid_level: int = 3,
    conv_tol: float = 1.0e-9,
    max_cycle: int = 100,
    verbose: int = 0,
    cube_points: tuple[int, int, int] = (80, 80, 80),
    cube_resolution: float | None = None,
    cube_margin: float = 5.0,
    molden_ignore_high_l: bool = True,
) -> FTDFTResult:
    """Run molecular finite-temperature KS-DFT and an FOD analysis.

    Args:
        mol: A built PySCF molecule.
        xc: Any XC expression accepted by PySCF.
        temperature: Fictitious electronic temperature in kelvin.  At zero,
            ordinary integer-occupation KS-DFT is performed.
        output_prefix: Prefix for ``.molden`` and ``.fod.cube`` files.  If
            omitted, no visualization files are written.
        unrestricted: Use UKS.  The default selects UKS for nonzero molecular
            spin and RKS otherwise.
        grid_level: PySCF numerical integration-grid level.
        conv_tol: SCF energy convergence tolerance in hartree.
        max_cycle: Maximum SCF cycles.
        verbose: PySCF output verbosity.
        cube_points: Grid dimensions used by :func:`pyscf.tools.cubegen.density`.
        cube_resolution: Optional cube spacing in bohr; when supplied it
            overrides ``cube_points`` as in PySCF.
        cube_margin: Empty space around the molecule in bohr.
        molden_ignore_high_l: Let PySCF omit unsupported angular momenta
            (l >= 5) from the Molden output.
    """
    if not isinstance(mol, gto.Mole) or not mol._built:
        raise TypeError("mol must be a built pyscf.gto.Mole")
    sigma = kelvin_to_hartree(temperature)
    if unrestricted is None:
        unrestricted = mol.spin != 0
    if not unrestricted and mol.spin != 0:
        raise ValueError("a molecule with nonzero spin requires unrestricted=True")
    if grid_level < 0:
        raise ValueError("grid_level must be non-negative")
    if max_cycle < 1:
        raise ValueError("max_cycle must be positive")
    if len(cube_points) != 3 or any(point < 2 for point in cube_points):
        raise ValueError("cube_points must contain three integers of at least 2")
    if cube_resolution is not None and cube_resolution <= 0.0:
        raise ValueError("cube_resolution must be positive")
    if cube_margin <= 0.0:
        raise ValueError("cube_margin must be positive")

    mean_field = dft.UKS(mol) if unrestricted else dft.RKS(mol)
    mean_field.xc = xc
    mean_field.grids.level = grid_level
    mean_field.conv_tol = conv_tol
    mean_field.max_cycle = max_cycle
    mean_field.verbose = verbose

    if sigma > 0.0:
        mean_field = scf.addons.smearing(
            mean_field,
            sigma=sigma,
            method="fermi",
            fix_spin=bool(unrestricted),
        )

    internal_energy = float(mean_field.kernel())
    if not mean_field.converged:
        raise RuntimeError(f"FT-DFT did not converge in {max_cycle} SCF cycles")

    entropy = float(mean_field.entropy) if sigma > 0.0 else 0.0
    free_energy = float(mean_field.e_free) if sigma > 0.0 else internal_energy
    zero_temperature_energy = (
        float(mean_field.e_zero) if sigma > 0.0 else internal_energy
    )
    mo_energy = np.asarray(mean_field.mo_energy, dtype=float).copy()
    mo_occ = np.asarray(mean_field.mo_occ, dtype=float).copy()
    weights = fod_weights(mo_occ, unrestricted=unrestricted)
    n_fod = float(np.sum(weights))
    fod_dm = fod_density_matrix(
        mean_field.mo_coeff, mean_field.mo_occ, unrestricted=unrestricted
    )

    # Since canonical MOs are overlap-orthonormal, Tr[D_FOD S] must equal the
    # sum of FOD weights.  This catches coefficient/occupation shape mistakes.
    ao_integral = float(np.einsum("ij,ji->", fod_dm, mean_field.get_ovlp()))
    if not math.isclose(ao_integral, n_fod, rel_tol=1.0e-8, abs_tol=1.0e-8):
        raise RuntimeError(
            f"inconsistent FOD density matrix: integral={ao_integral}, N_FOD={n_fod}"
        )

    molden_file: Path | None = None
    fod_cube_file: Path | None = None
    cube_integral: float | None = None
    if output_prefix is not None:
        molden_file, fod_cube_file = _output_paths(output_prefix)
        molden.from_scf(
            mean_field, str(molden_file), ignore_h=molden_ignore_high_l
        )
        density = cubegen.density(
            mol,
            str(fod_cube_file),
            fod_dm,
            nx=cube_points[0],
            ny=cube_points[1],
            nz=cube_points[2],
            resolution=cube_resolution,
            margin=cube_margin,
        )
        cube_integral = _cube_integral(
            mol, density, cube_points, cube_resolution, cube_margin
        )
        allowed_error = max(5.0e-3, 0.01 * n_fod)
        if abs(cube_integral - n_fod) > allowed_error:
            warnings.warn(
                "The finite cube integrates to "
                f"{cube_integral:.8f}, while N_FOD is {n_fod:.8f}. "
                "Increase cube_margin or refine the cube grid.",
                RuntimeWarning,
                stacklevel=2,
            )

    return FTDFTResult(
        temperature_kelvin=float(temperature),
        sigma_hartree=sigma,
        internal_energy=internal_energy,
        free_energy=free_energy,
        zero_temperature_energy=zero_temperature_energy,
        entropy_kb=entropy,
        n_fod=n_fod,
        cube_integral=cube_integral,
        mo_energy=mo_energy,
        mo_occ=mo_occ,
        unrestricted=unrestricted,
        molden_file=molden_file,
        fod_cube_file=fod_cube_file,
        mean_field=mean_field,
    )
