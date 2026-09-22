# N2 FOD comparison: ORCA and PySCF

This comparison checks whether the PySCF implementation reproduces the
practical `N_FOD` trend from ORCA for N2 bond stretching. It is deliberately a
default-settings comparison rather than an attempt to make every numerical
integration and SCF option identical.

## Common model chemistry

- molecule: neutral singlet N2;
- bond lengths: 1.10, 1.40, 1.80, 2.20, 2.60, and 3.00 angstrom;
- functional: TPSS;
- orbital basis: def2-TZVP;
- electronic temperature: 5000 K;
- restricted finite-temperature KS-DFT.

ORCA 6.1.0 was run using only `! FOD`, which selects its documented FOD
defaults. In this installation ORCA also selected its default RI-J treatment
with the def2/J auxiliary basis. PySCF 2.13.0 was run with
`run_ftdft(..., xc="tpss", temperature=5000.0)` and otherwise used the package
defaults. These program-specific numerical defaults were not forced to match.

## Results

| R(N-N) / A | ORCA N_FOD | PySCF N_FOD | PySCF - ORCA |
| ---: | ---: | ---: | ---: |
| 1.10 | 0.000260 | 0.000260 | +0.000000 |
| 1.40 | 0.020115 | 0.020124 | +0.000009 |
| 1.80 | 0.564187 | 0.564166 | -0.000021 |
| 2.20 | 2.010567 | 2.009684 | -0.000883 |
| 2.60 | 3.388195 | 3.388107 | -0.000088 |
| 3.00 | 4.464895 | 4.465076 | +0.000181 |

Over these six geometries:

- mean absolute difference: `0.000197`;
- root-mean-square difference: `0.000370`;
- maximum absolute difference: `0.000883` at 2.20 angstrom.

The two programs give the same physical picture. `N_FOD` is essentially zero
near the equilibrium bond length and rises strongly as the N2 triple bond is
stretched. The small residual differences are consistent with leaving each
program's SCF, density-fitting, and numerical-grid defaults in place. ORCA
prints `N_FOD` to six decimal places, which also limits the precision of the
reported comparison.

The machine-readable values are in
[`n2_fod_comparison.csv`](n2_fod_comparison.csv).

## Reproduce

After installing this repository, run:

```bash
python comparisons/n2_fod/run_comparison.py \
  --orca /full/path/to/orca \
  --output-dir n2_fod_runs
```

The script generates a minimal ORCA input for every distance, runs ORCA and
PySCF, parses the ORCA `N_FOD` value, and writes a new comparison CSV. A
representative ORCA input is provided in
[`orca_default_example.inp`](orca_default_example.inp).

