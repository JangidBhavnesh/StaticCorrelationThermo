from pathlib import Path

import numpy as np
import pytest
from pyscf import gto
from pyscf.tools import molden

from static_correlation_thermo import (
    BOLTZMANN_HARTREE_PER_KELVIN,
    fod_density_matrix,
    fod_weights,
    kelvin_to_hartree,
    run_ftdft,
)


def test_temperature_conversion_and_validation():
    assert kelvin_to_hartree(1000.0) == pytest.approx(
        1000.0 * BOLTZMANN_HARTREE_PER_KELVIN
    )
    with pytest.raises(ValueError):
        kelvin_to_hartree(-1.0)


def test_restricted_fod_weights_and_density_matrix():
    occupations = np.array([2.0, 1.5, 0.5, 0.0])
    expected = np.array([0.0, 0.5, 0.5, 0.0])
    assert np.allclose(fod_weights(occupations), expected)
    assert np.allclose(fod_density_matrix(np.eye(4), occupations), np.diag(expected))


def test_unrestricted_fod_weights_and_density_matrix():
    occupations = np.array([[1.0, 0.25], [0.75, 0.0]])
    coefficients = np.array([np.eye(2), np.eye(2)])
    expected_weights = np.array([[0.0, 0.25], [0.25, 0.0]])
    assert np.allclose(
        fod_weights(occupations, unrestricted=True), expected_weights
    )
    assert np.allclose(
        fod_density_matrix(coefficients, occupations, unrestricted=True),
        np.diag([0.25, 0.25]),
    )


def test_stretched_h2_writes_consistent_outputs(tmp_path: Path):
    molecule = gto.M(
        atom="H 0 0 -1.5; H 0 0 1.5",
        basis="sto-3g",
        unit="Angstrom",
        verbose=0,
    )
    result = run_ftdft(
        molecule,
        xc="lda,vwn",
        temperature=5000.0,
        output_prefix=tmp_path / "h2",
        grid_level=0,
        conv_tol=1.0e-8,
        cube_points=(32, 32, 32),
        cube_margin=5.0,
    )

    assert result.molden_file == tmp_path / "h2.molden"
    assert result.fod_cube_file == tmp_path / "h2.fod.cube"
    assert result.molden_file.is_file()
    assert result.fod_cube_file.is_file()
    assert result.entropy_kb > 0.0
    assert result.n_fod > 0.1
    assert result.free_energy == pytest.approx(
        result.internal_energy - result.sigma_hartree * result.entropy_kb
    )
    assert result.cube_integral == pytest.approx(result.n_fod, abs=5.0e-3)

    _, _, _, molden_occ, _, _ = molden.load(str(result.molden_file))
    assert np.asarray(molden_occ) == pytest.approx(result.mo_occ)


def test_open_shell_uses_uks_and_preserves_spin_populations():
    molecule = gto.M(
        atom="Li 0 0 0", basis="sto-3g", spin=1, verbose=0
    )
    result = run_ftdft(
        molecule,
        xc="lda,vwn",
        temperature=5000.0,
        grid_level=0,
        conv_tol=1.0e-8,
    )

    assert result.unrestricted
    assert np.sum(result.mo_occ, axis=1) == pytest.approx(molecule.nelec)
    assert result.n_fod == pytest.approx(
        np.sum(fod_weights(result.mo_occ, unrestricted=True))
    )
