# Parametric Beta-Barrel / Solenoid Backbone Generation

This module provides parametric backbone generation for beta-barrels, solenoids (beta-helices), and the strand components of TIM barrels. It extends the `core::conformation::parametric` framework with a new `BarrelParametrizationCalculator` and associated movers, following the same architecture as `protocols/helical_bundle/`.

## Geometric Model

Each strand is placed on a cylindrical (or helical) surface using the generalized Crick equations from `numeric::crick_equations`. The per-atom geometry of a beta strand comes from `database/protocol_data/crick_parameters/beta_strand.crick_params`. The barrel/solenoid shape is controlled by:

- **`r0`** — radius of the cylinder (distance from barrel axis to strand backbone)
- **`omega0`** — superhelical twist per residue. **This is the key parameter distinguishing barrels from solenoids**: `omega0 = 0` produces a closed barrel (strands parallel to the axis); `omega0 != 0` produces a solenoid where strands spiral around the axis.
- **`n_strands`** / **`shear_number`** — global barrel topology parameters. The shear number `S` controls the axial stagger between adjacent strands: strand `j` is offset by `(j-1) * S * z1 / n` along the barrel axis.
- **`antiparallel`** — if true, every other strand is inverted (as in OMP barrels). If false, all strands run in the same direction (as in TIM barrels).
- **`epsilon`** — lateral squash for non-circular cross-sections.
- **`delta_omega0`** / **`delta_z0`** — per-strand azimuthal and axial offsets (computed automatically from `n_strands` and `shear_number`, but can be overridden).

## Movers

### MakeBarrel

Builds a complete parametric barrel or solenoid from scratch.

```xml
<!-- 8-strand all-parallel barrel (TIM barrel inner sheet) -->
<MakeBarrel name="make_barrel" n_strands="8" shear_number="8"
    r0="8.0" antiparallel="false" use_degrees="false"
    crick_params_file="beta_strand" residue_name="ALA"
    strand_length="10" reset="true" />

<!-- 16-strand antiparallel barrel (OMP-like) -->
<MakeBarrel name="make_omp" n_strands="16" shear_number="20"
    r0="12.0" antiparallel="true"
    crick_params_file="beta_strand" residue_name="ALA"
    strand_length="8" />

<!-- Solenoid / beta-helix (nonzero omega0) -->
<MakeBarrel name="make_solenoid" n_strands="3" shear_number="0"
    r0="7.0" omega0="0.05" antiparallel="false"
    crick_params_file="beta_strand" residue_name="ALA"
    strand_length="15" />
```

Per-strand overrides can be specified with `<Strand>` sub-tags:

```xml
<MakeBarrel name="make" n_strands="8" shear_number="8" r0="8.0"
    antiparallel="false" crick_params_file="beta_strand"
    residue_name="ALA" strand_length="10">
  <Strand strand_length="12" />  <!-- strand 1: longer -->
  <Strand strand_length="8" />   <!-- strand 2: shorter -->
  <!-- strands 3-8 use the default strand_length="10" -->
</MakeBarrel>
```

### PerturbBarrel

Perturbs barrel parameters during Monte Carlo sampling. Supports Gaussian and uniform perturbation types.

```xml
<PerturbBarrel name="perturb" default_perturbation_type="gaussian"
    use_degrees="false" r0_perturbation="0.3">
  <!-- Per-strand overrides -->
  <Strand strand_index="1" delta_omega0_perturbation="0.05" />
  <Strand strand_index="3" delta_z0_perturbation="0.2" />
</PerturbBarrel>
```

## Configurations for Specific Structure Types

### Closed beta-barrels (OMP family)

Outer membrane protein barrels are antiparallel, closed cylinders with 8-22 strands.

- `omega0 = 0` (strands parallel to barrel axis)
- `antiparallel = true`
- Typical `r0`: 7-15 A depending on strand count
- Shear numbers are typically `S = n` or `S = n + 2`
- Strand lengths 6-25 residues (variable per strand in real OMPs)

### TIM barrels (strand component)

TIM barrels have 8 parallel beta-strands forming an inner barrel, with alpha-helices connecting them on the outside.

- `omega0 = 0`, `antiparallel = false`, `n_strands = 8`, `shear_number = 8`
- Typical `r0`: 7-9 A
- For the alpha-helical components, use the `MakeBundle` / `MakeBundleHelix` movers from `protocols/helical_bundle/` with a larger radius

### Solenoids / beta-helices

Beta-helices are open, non-closing structures where strands spiral around the axis. Each "rung" typically has 2-3 short strands at different azimuthal positions.

- `omega0 != 0` (the defining difference from closed barrels)
- `antiparallel = false` (typically parallel)
- Smaller `n_strands` (2-3 per rung) with longer `strand_length`
- `shear_number = 0` (stagger is provided by the non-zero omega0)
- Adjust `omega0` to control the pitch of the solenoid

### Collagen triple helices

Collagen does **not** use this module. It is already supported by `MakeBundle` with `crick_params_file="collagen"` and 3-fold symmetry. See the integration test `make_and_perturb_bundle_multirepeat` for an example.

## Class Hierarchy

```
core::conformation::parametric::ParametrizationCalculator  (base, in core/)
  └── BarrelParametrizationCalculator                      (this module)

core::conformation::parametric::Parameters                 (base, in core/)
  └── BarrelParameters                                     (per-strand container)

core::conformation::parametric::ParametersSet              (base, in core/)
  └── BarrelParametersSet                                  (whole-barrel container)
```

## Relationship to helical_bundle

This module reuses the coordinate generation infrastructure from `protocols/helical_bundle/util.hh` (`generate_atom_positions`, `place_atom_positions`, etc.) and the Crick equation math from `numeric/crick_equations/`. The `.crick_params` file format is shared. The architectural pattern (Calculator + Parameters + ParametersSet + Make/Perturb movers) is identical.
