# H2O HF/6-31G(d) reference wavefunctions for QTAIM and density tests

The same water molecule at one fixed geometry, computed three ways, with the
wavefunction written in every format Avogadro (or its planned readers) needs,
plus point values of the electron density and its derivatives from critic2.

## Geometry

`geometry.xyz` (Angstrom). It is the B3LYP/def2-TZVP optimized structure from
ORCA 6.1.0 (`! B3LYP/G def2-TZVP NoRI VeryTightSCF TightOpt Freq`), kept fixed
for every calculation here. The Gaussian runs use `nosymm`, so the molecule
is not reoriented and all files share this frame.

`points.xyz` holds 20 comparison points (Angstrom, same frame): 2 bond
midpoints, 14 valence-region points and 4 tail points 2.5 to 5 Bohr from the
nearest nucleus.

## Calculations

All are restricted Hartree-Fock, closed shell, charge 0, singlet, with the
6-31G(d) basis. 6D and 5D are different basis sets, not different encodings
of one basis, so their densities differ; compare a file only with data from the
same basis.

| Folder | Program | d functions | Files |
|--------|---------|-------------|-------|
| `gaussian-hf-6-31gd-6d/` | Gaussian 16 Rev. C.01 | Cartesian (6D) | `h2o.fchk.gz`, `h2o.wfn`, `h2o.wfx`, `h2o.gjf` |
| `gaussian-hf-6-31gd-5d/` | Gaussian 16 Rev. C.01 | spherical (5D) | `h2o.fchk.gz`, `h2o.wfn`, `h2o.wfx`, `h2o.gjf` |
| `orca-hf-6-31gd-5d/` | ORCA 6.1.0 | spherical (5D) | `h2o.molden`, `h2o.wfn`, `h2o.wfx`, `h2o.inp` |

- Gaussian: `HF/6-31G(d) 6D 10F` or `5D 7F` as appropriate, `nosymm
  scf=tight int=ultrafine`. The `.wfx` comes from the first job step and the
  `.wfn` from a second step reading the checkpoint (see `h2o.gjf`). The `.fchk`
  is gzip compressed (as in `data/integration/`).
- ORCA: `! HF 6-31G(d) AIM NoRI VeryTightSCF` (see `h2o.inp`). The Molden file
  was written by `orca_2mkl`; the `.wfn` and `.wfx` by the `AIM` keyword.
- Program output logs are deliberately not included.

## reference-points.json

For every wavefunction file above, critic2 1.4.306 (using the file itself):
density `rho`, `gradient`, `hessian`, `laplacian`, and the kinetic energy
densities `G` and `K` at each of the 20 points. Point positions are given in
both Angstrom (verbatim from `points.xyz`) and Bohr (as critic2 used them,
converted with 0.529177210903 A/Bohr). Units are atomic units; Bohr for
lengths. Values are copied exactly as critic2 printed them. Provenance for
each entry (the file critic2 read and its summary) is recorded in the JSON.

Regenerate with
`python3 make_reference_points.py <path to inputs/h2o> reference-points.json`;
the script needs the critic2 `summary.json` files (schema
`avogadro-qtaim-summary/1`), which are not stored here.

## License

These files are distributed with the `avogadrodata` repository. That
repository has no separate LICENSE file; its README states that it is
distributed under the 3-clause BSD license (BSD-3-Clause), the same as
Avogadro. The input files and `reference-points.json` are original to this
project.
