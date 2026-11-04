# H2O HF/6-31G(d) reference wavefunctions for QTAIM and density tests

The same water molecule at one fixed geometry, computed with Gaussian, ORCA,
Psi4, PySCF and xtb, with the wavefunction written in every format Avogadro
(or its planned readers) needs, plus point values of the electron density and
its derivatives from critic2.

## Geometry

`geometry.xyz` (Angstrom). It is the B3LYP/def2-TZVP optimized structure from
ORCA 6.1.0 (`! B3LYP/G def2-TZVP NoRI VeryTightSCF TightOpt Freq`), kept fixed
for every calculation here. The Gaussian runs use `nosymm`, so the molecule
is not reoriented and all files share this frame.

`points.xyz` holds 20 comparison points (Angstrom, same frame): 2 bond
midpoints, 14 valence-region points and 4 tail points 2.5 to 5 Bohr from the
nearest nucleus.

## Calculations

All but xtb are restricted Hartree-Fock, closed shell, charge 0, singlet, with
the 6-31G(d) basis. 6D and 5D are different basis sets, not different encodings
of one basis, so their densities differ; compare a file only with data from the
same basis.

| Folder | Program | d functions | Files |
|--------|---------|-------------|-------|
| `gaussian-hf-6-31gd-6d/` | Gaussian 16 Rev. C.01 | Cartesian (6D) | `h2o.fchk.gz`, `h2o.wfn`, `h2o.wfx`, `h2o.gjf` |
| `gaussian-hf-6-31gd-5d/` | Gaussian 16 Rev. C.01 | spherical (5D) | `h2o.fchk.gz`, `h2o.wfn`, `h2o.wfx`, `h2o.gjf` |
| `orca-hf-6-31gd-5d/` | ORCA 6.1.0 | spherical (5D) | `h2o.molden`, `h2o.wfn`, `h2o.wfx`, `h2o.inp` |
| `psi4-hf-6-31gd-6d/` | Psi4 1.11 | Cartesian (6D) | `h2o.molden`, `h2o.fchk.gz`, `h2o.in` |
| `psi4-hf-6-31gd-5d/` | Psi4 1.11 | spherical (5D) | `h2o.molden`, `h2o.fchk.gz`, `h2o.in` |
| `pyscf-hf-6-31gd-6d/` | PySCF 2.14.0 | Cartesian (6D) | `h2o.molden`, `h2o.py` |
| `pyscf-hf-6-31gd-5d/` | PySCF 2.14.0 | spherical (5D) | `h2o.molden`, `h2o.py` |
| `xtb-gfn2/` | xtb 6.7.1 | Cartesian d (minimal valence basis) | `h2o.molden`, `run.sh` |

- Gaussian: `HF/6-31G(d) 6D 10F` or `5D 7F` as appropriate, `nosymm
  scf=tight int=ultrafine`. The `.wfx` comes from the first job step and the
  `.wfn` from a second step reading the checkpoint (see `h2o.gjf`). The `.fchk`
  is gzip compressed (as in `data/integration/`).
- ORCA: `! HF 6-31G(d) AIM NoRI VeryTightSCF` (see `h2o.inp`). The Molden file
  was written by `orca_2mkl`; the `.wfn` and `.wfx` by the `AIM` keyword.
- Psi4: restricted HF/6-31G(d), see `h2o.in` (its header comment says it was
  written for Psi4 1.3.2; the run that produced these files was Psi4 1.11).
  The Molden and fchk files come from Psi4 itself.
- PySCF: restricted HF/6-31G(d) with `cart=True` (6D) or `cart=False` (5D),
  see `h2o.py`; the Molden file is written with `pyscf.tools.molden`.
- xtb: GFN2-xTB single point at the same geometry (`run.sh`,
  `xtb --gfn 2 --molden`). Its density is valence only (a minimal valence
  basis, no core), so it is not comparable to the Hartree-Fock values; the
  critic2 numbers were computed from the Molden file only.
- Program output logs are deliberately not included.

Conventions of the Molden files (what `QuantumIO::MoldenFile` has to cope
with): Psi4 and PySCF write the coefficients of normalized primitives, the
convention Avogadro assumes. xtb includes the normalization of
x^l exp(-alpha r^2) in every coefficient, and ORCA (`orca_2mkl`) does the same
with a different angular factor for d, f and g, and also negates some of the
spherical f and g components; see `GaussianSet::PrimitiveNormalization`.

## reference-points.json

For every wavefunction file above, critic2 1.4.306 (using the file itself, in
the `critic2_dialect` recorded with each Molden entry: `orca` for ORCA,
`pyscf` for PySCF and also for Psi4 1.11, whose Molden files critic2's default
and `psi4` dialects read about 1.8e-2 off for Cartesian d, and `xtb` for xtb):
density `rho`, `gradient`, `hessian`, `laplacian`, and the kinetic energy
densities `G` and `K` at each of the 20 points. Point positions are given in
both Angstrom (verbatim from `points.xyz`) and Bohr (as critic2 used them,
converted with 0.529177210903 A/Bohr). Units are atomic units; Bohr for
lengths. Values are copied exactly as critic2 printed them. Provenance for
each entry (the file critic2 read and its summary) is recorded in the JSON.

The xtb entry has no QTAIM run (the density is valence only): its values are
the point values (`rho`, `gradient`, `hessian`, `laplacian`, `G`; no `K`)
critic2 printed in `%% point` blocks of its output, with `load ... xtb`, which
are not stored here.

Regenerate with `python3 make_reference_points.py <path to inputs>`, which
also writes the files of the neighbouring `h2o-hf-ccpvtz/` and
`n2-b3lyp-ccpvqz/` folders (and, given a second argument, `../molden/
SF4-reference-points.json`); the script needs the critic2 `summary.json` files
(schema `avogadro-qtaim-summary/1`), which are not stored here.

## License

These files are distributed with the `avogadrodata` repository. That
repository has no separate LICENSE file; its README states that it is
distributed under the 3-clause BSD license (BSD-3-Clause), the same as
Avogadro. The input files and `reference-points.json` are original to this
project.

## SF4 (`data/molden/SF4-reference-points.json`)

`data/molden/SF4.molden` is an xtb Molden file with Cartesian d functions on
sulfur. `data/molden/SF4-reference-points.json` holds critic2 1.4.306 point
values (`load SF4.molden xtb`, `rho`, `gradient`, `hessian`, `laplacian`, `G`)
at 12 points around the sulfur atom, in Angstrom in the file's own frame.
They are written by `make_reference_points.py` given the critic2 point output.
