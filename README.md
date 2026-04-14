# DFT from Scratch — A Python Tutorial

A self-contained Jupyter notebook that implements Kohn-Sham Density Functional Theory
from first principles using only NumPy and SciPy. Every step — basis set construction,
numerical integration, the Hartree and exchange-correlation potentials, and the SCF loop —
is written explicitly so the physics remains visible in the code.

**Author:** Marco Fronzi

---

## Overview

Most DFT tutorials either stop at the equations or jump straight to a production code
(VASP, Quantum ESPRESSO, CP2K) where the implementation details are hidden. This
notebook occupies the middle ground: the code is intentionally simple and slow, trading
computational efficiency for transparency. After working through it you will understand
what a DFT engine actually computes at each step, and why.

The notebook is structured as a linear narrative. Markdown cells derive or explain each
equation; the code cell immediately following implements it. You can run each section
independently once the setup cells have been executed.

---

## What You Will Learn

- The Hohenberg-Kohn theorems and why the ground-state energy is a functional of the density
- The Kohn-Sham mapping: replacing the interacting many-body problem with a set of
  single-particle equations in an effective potential
- How atomic Slater-type orbitals (STOs) are constructed from Gaussian primitives and
  spherical harmonics
- How a double-zeta (DZ) basis set is designed so that the occupied orbital space has
  genuine degrees of freedom (the `n_virt ≥ n_occ` criterion)
- Numerical evaluation of the one-electron matrices: kinetic energy **T**, nuclear
  attraction **V**, and overlap **S**, using finite-difference Laplacians and real-space
  grid integration
- The Hartree potential `v_H(r)` as a pairwise Coulomb sum over the density and its
  contribution to the Fock matrix
- The LDA Slater exchange functional: energy `E_xc` and potential `v_xc = δE_xc/δρ`
- Construction of the Kohn-Sham Fock matrix: `F = H_core + J[ρ] + V_xc[ρ]`
- Canonical orthogonalisation of the overlap matrix to handle near-linear dependence
- The Roothaan-Hall generalised eigenvalue problem `FC = SCε` and its solution
- The SCF cycle: density mixing, energy convergence criterion, and electron-count
  conservation `Tr(PS) = N_el`
- Post-SCF analysis: orbital energy spectrum, HOMO/LUMO identification, and charge
  density visualisation

---

## Notebook Structure

| Section | Cells | Content |
|---------|-------|---------|
| **Theory** | 00–06 | Hohenberg-Kohn theorems, KS equations, variational derivation, LCAO expansion, STO basis function form |
| **Setup** | 07–10 | Imports; `ELEMENTS` table; SZ and DZ basis set parameters; core classes (`Atom`, `GaussianBasisFunction`, `MolecularOrbital`); molecule geometry and EHT initial guess |
| **Pre-SCF visualisation** | 11–16 | 2-D and 3-D plots of basis functions and initial molecular orbitals |
| **XC functionals** | 17–19 | LDA Slater exchange and simplified PBE-like GGA; energy and potential |
| **Matrix elements** | 20–25 | Overlap **S**, kinetic **T**, nuclear attraction **V**, Hartree potential and matrix, LDA XC matrix; initial density matrix |
| **Fock matrix & solver** | 26–28 | `compute_fock_matrix`, `solve_roothaan_hall` (canonical orthogonalisation), `compute_density_matrix`, density mixing |
| **SCF loop** | 29–31 | `compute_total_energy`, `scf_cycle` with energy convergence and `Tr(PS)` diagnostic; run cell |
| **Post-SCF analysis** | 32–34 | Orbital energy table (Ha and eV), HOMO-LUMO gap, converged MO coefficients, charge density; energy level diagram; 2-D HOMO/LUMO/density plots |

---

## Key Implementation Decisions

### Basis set — why DZ instead of minimal SZ

The notebook ships with two basis sets:

| Basis | Functions per atom | H₂O total | n_occ | n_virt |
|-------|--------------------|-----------|-------|--------|
| `SZ` (minimal) | 1 per shell | 6 | 5 | **1** |
| `DZ` (double-zeta, **default**) | 2 per shell + core | 13 | 5 | **8** |

With `SZ`, the canonical orthogonalisation discards the single near-dependent vector,
leaving exactly five usable vectors for five occupied orbitals. The SCF has no degrees
of freedom and converges trivially in one step from any starting density. The `DZ`
basis adds an inner and outer exponent for every shell (plus an explicit 1s core for
second-row atoms), giving `n_virt ≥ n_occ` for all molecules in the tutorial and
genuine multi-step SCF iteration.

### Convergence criterion — energy, not `‖ΔP‖`

In a non-orthogonal AO basis, `Tr(P) ≫ N_el` while `Tr(PS) = N_el`. The Frobenius
norm `‖P_new − P_old‖_F` is dominated by the basis metric, not the electron density,
and is not a reliable convergence measure. The SCF loop converges on `|E_new − E_old|`,
which is the physically meaningful criterion used in production codes.

### Hartree potential — pairwise Coulomb sum

The Hartree potential is built as a direct pairwise sum:

```
v_H(r_i) = Σ_j  ρ(r_j) / |r_i − r_j|  · dV
```

with diagonal regularised by the Wigner-Seitz cell radius `r_cut = dV^(1/3)` to avoid
self-interaction divergence. The matrix element `J_ij = ∫ φ_i(r) v_H(r) φ_j(r) dr`
is then evaluated by numerical quadrature. This is `O(N_grid²)` and therefore slow for
large grids, but pedagogically correct — every term maps directly to a formula in the
theory section.

### Total energy — no double-counting

```
E_tot = Tr[P H_core] + E_H + E_xc + E_nn
```

Using `E = ½ Tr[P(H_core + F)]` would double-count the Hartree and XC contributions
because both already appear in F. The correct decomposition is used throughout.

---

## Supported Molecules

Pre-configured geometries are included for:

| Molecule | Formula | Electrons | DZ basis functions |
|----------|---------|-----------|-------------------|
| Molecular hydrogen | H₂ | 2 | 4 |
| Water | H₂O | 10 | 13 |
| Ammonia | NH₃ | 10 | 15 |
| Molecular nitrogen | N₂ | 14 | 18 |
| Carbon monoxide | CO | 14 | 18 |
| Methane | CH₄ | 10 | 17 |

To switch molecule, uncomment the relevant geometry block in the setup cell. Any
closed-shell molecule built from H, He, Li–Ne, Na, Si, S, Cl can be added by defining
its Cartesian coordinates.

---

## Requirements

```
python >= 3.9
numpy
scipy
matplotlib
```

No external chemistry libraries are required. Install with:

```bash
pip install numpy scipy matplotlib
```

---

## Running the Notebook

```bash
jupyter notebook DFT_tutorial.ipynb
```

Run cells top-to-bottom. The SCF cell (the `scf_cycle` call) is the only slow step:
with the default 10-point-per-axis grid (1000 grid points) and H₂O/DZ (13 basis
functions) it takes roughly 1–3 minutes on a modern laptop. Reduce `n_points` to 8
(512 grid points) for faster iteration during exploration.

### Grid size vs accuracy

| `n_points` | Grid points | H₂O runtime (approx) | Use case |
|------------|-------------|----------------------|----------|
| 6 | 216 | ~15 s | Quick exploration |
| 8 | 512 | ~45 s | Development |
| 10 | 1 000 | ~2 min | Tutorial default |
| 15 | 3 375 | ~15 min | Better accuracy |

---

## Function Reference

### Basis set and geometry

| Function | Description |
|----------|-------------|
| `define_basis_set_params()` | Returns SZ and DZ exponent tables for H–Cl |
| `generate_atomic_orbitals(geom, params, name)` | Builds list of `GaussianBasisFunction` for a molecule |
| `initialize_molecular_orbitals(basis, n_el, method)` | EHT/AO/SAD initial MO coefficient matrix |
| `create_grid_points(x, y, z, n)` | Flat `(N³, 3)` array of Cartesian grid points |

### Integral engine

| Function | Description |
|----------|-------------|
| `precompute_basis_values(basis, grid)` | `(n_basis, n_grid)` matrix of φ_μ(r_p) |
| `kinetic_matrix(basis, grid)` | **T**: `−½∫φ_i ∇²φ_j dr` via numerical Laplacian |
| `nuclear_attraction_matrix(basis, atoms, grid)` | **V**: `∫φ_i V_nuc φ_j dr` |
| `overlap_matrix(basis, grid)` | **S**: `∫φ_i φ_j dr` |
| `compute_electron_density_on_grid(P, bfv)` | ρ(r) from density matrix |
| `compute_hartree_potential_on_grid(rho, grid)` | v_H(r) pairwise Coulomb sum |
| `hartree_matrix_contribution(vH, bfv, grid)` | **J**: `∫φ_i v_H φ_j dr` |
| `lda_xc_potential_on_grid(rho)` | v_xc(r) = −(4/3) C_X ρ^(1/3) |
| `xc_matrix_contribution(vxc, bfv, grid)` | **V_xc**: `∫φ_i v_xc φ_j dr` |

### SCF loop

| Function | Description |
|----------|-------------|
| `compute_fock_matrix(basis, P, atoms, grid, T, V, bfv)` | F = H_core + J + V_xc |
| `solve_roothaan_hall(F, S)` | FC = SCε via canonical orthogonalisation |
| `compute_density_matrix(C, n_el)` | P = 2 C_occ C_occ^T |
| `mix_density_matrices(P_old, P_new, alpha)` | Linear damping P ← (1−α)P_old + α P_new |
| `compute_total_energy(P, T, V, rho, vH, grid, atoms)` | E = Tr[PH_core] + E_H + E_xc + E_nn |
| `scf_cycle(basis, atoms, grid, C0, n_el, ...)` | Full SCF loop with convergence reporting |

### Post-SCF visualisation

| Function | Description |
|----------|-------------|
| `plot_homo_lumo_orbital_diagram(energies, n_el)` | Vertical energy level diagram |
| `plot_homo_lumo_density(C, P, basis, atoms, ...)` | Three-panel 2-D slice: HOMO, LUMO, ρ(r) |
| `plot_atomic_orbitals_2d(basis, plane, ...)` | Individual basis function plots |
| `plot_molecular_orbitals_2d(mos, plane, ...)` | Initial MO wavefunction plots |

---

## Limitations and Pedagogical Caveats

This is a tutorial implementation, not a production code. Known limitations:

- **Numerical integration only.** All integrals are evaluated on a real-space grid.
  Analytic Gaussian integral formulas (Boys functions, McMurchie-Davidson recursion)
  are not used. Accuracy improves with grid density but convergence is slow.

- **No periodic boundary conditions.** Molecules only. Solid-state DFT requires Bloch
  theorem, k-point sampling, and plane-wave or PAW bases.

- **LDA exchange only.** Correlation is not included. A full LDA functional would add
  the Vosko-Wilk-Nusair (VWN) or Perdew-Wang (PW92) correlation. GGA (PBE) is
  implemented but not used in the SCF by default.

- **Closed-shell only.** The density matrix uses double occupancy throughout. Open-shell
  systems require spin-unrestricted (UKS) or spin-restricted open-shell (ROKS)
  formulations.

- **Illustrative exponents.** The STO exponents in the DZ basis are calibrated to
  reproduce qualitative orbital shapes, not optimised for energy accuracy. Results are
  not comparable to production 6-31G or cc-pVDZ calculations.

- **O(N_grid²) Hartree term.** The pairwise Coulomb sum scales quadratically with the
  number of grid points. Production codes use the Poisson equation (FFT-based) for O(N
  log N) scaling.

---

## Background Reading

- W. Kohn and L. J. Sham, *Self-Consistent Equations Including Exchange and Correlation
  Effects*, Phys. Rev. **140**, A1133 (1965)
- P. Hohenberg and W. Kohn, *Inhomogeneous Electron Gas*, Phys. Rev. **136**, B864 (1964)
- R. G. Parr and W. Yang, *Density-Functional Theory of Atoms and Molecules*, Oxford (1989)
- C. J. Cramer, *Essentials of Computational Chemistry*, 2nd ed., Wiley (2004) — Chapter 8
- K. Burke, *The ABC of DFT*, lecture notes (2007) —
  [online](https://dft.uci.edu/doc/g1.pdf)

