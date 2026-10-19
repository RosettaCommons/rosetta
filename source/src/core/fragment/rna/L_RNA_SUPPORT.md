# Mixed D-RNA / L-RNA modeling (FARFAR / rna_denovo)

This note covers what this fork changes so that D-RNA and L-RNA (mirror-image RNA) can be modeled
in the same pose, and what is still D-only.

## Short version

* Write L-RNA residues as `1`, `2`, `3`, `4` (L-a, L-c, L-g, L-u) in `-sequence` or the FASTA file.
  `a c g u` are still D-RNA. Example for a D strand and an L strand:

  ```
  >D chain A:1-8
  ggcgaagc
  >L chain B:1-8
  33232213
  ```

* No custom params files are needed. Rosetta generates the L-RNA residue types (`LRAD`, `LRCY`,
  `LRGU`, `LURA`, PDB names ` 0A`, ` 0C`, ` 0G`, ` 0U`) from the D types with the `L_RNA` patch
  (`patches/nucleic/rna/L_RNA.txt`). Their geometry is an exact mirror of the D types.
* Everything below turns on automatically when the pose contains L-RNA. In a pure D-RNA pose
  every changed code path reduces to the original arithmetic, so results are unchanged.
* Override: `-rna:denovo:invert_torsions_for_l_rna false` disables fragment and jump mirroring
  (the old behavior). Only use it if `-vall_torsions` points at a library that is already mirrored
  and the pose is pure L-RNA. Scoring stays chirality-aware either way.
* The first time an L-RNA window receives a fragment, the `protocols.rna.denovo.fragments.FullAtomRNA_Fragments`
  tracer prints one line:
  `L-RNA detected: inserting mirror-image (sign-inverted) vall fragments into L-RNA residues ...`

## How mirroring works

Reflecting a structure keeps every bond length and angle and negates every dihedral. The vall,
the jump library, the ideal-geometry PDBs and every knowledge-based potential in Rosetta come from
D-RNA. So an L-RNA residue is handled as the mirror image of a D-RNA residue:

| Quantity | L-RNA handling |
|---|---|
| Backbone torsions α–ζ, χ, ν2, ν1, 2'-OH | negated (`rna_chirality_factor()` = −1) |
| Non-main-chain sugar coordinates (C3'/C4'/C5' stub frame) | local `z` negated (`mirror_non_main_chain_sugar_coords()`) |
| Base-pair jumps (`RNA_JumpMover`) | `M → F·M·F`, `t → F·t`, `F = diag(1,1,−1)` |
| Base coordinate frame used for scoring | `z` axis negated (`get_rna_base_coordinate_system_chirality_normalized()`) |
| Ideal-geometry reference residues (`RNA_IdealCoord`) | mirrored copies (x → −x) of the D PDBs |

The stub-based rules follow from how Rosetta builds frames: the third axis is always a cross
product, so reflecting a residue keeps the local x and y coordinates and negates z.

## What changed, by file

### Detection and helpers
* `core/chemical/rna/util.{hh,cc}`: `rna_chirality_factor()`,
  `get_rna_base_coordinate_system_chirality_normalized()` and its `_stub_` variant. `get_residue_pucker_state()`
  and `get_residue_base_state()` (syn/anti) now classify the mirrored torsion.
* Chirality comes from `ResidueType::is_l_rna()` (the `L_RNA` property).

### Fragment insertion (vall)
* `core/fragment/rna/FullAtomRNA_Fragments.{hh,cc}`
  * `insert_fragment()` negates every torsion and mirrors the sugar coordinates when the window is
    L-RNA. All callers get this: `apply_random_fragment`, `apply_best_fragment`, and the symmetry copies.
  * `window_has_uniform_rna_chirality()`: a window that mixes D and L residues is never filled.
    `insert_fragment()` skips it with a one-time warning.
  * Vall lookup maps L residues to the D sequence (`1` → `a`, ...), so homology exclusion and
    sequence matching work unchanged.
* `protocols/rna/denovo/movers/RNA_FragmentMover.cc`: mixed-chirality windows are removed from the
  insertion map, so they are never picked and no Monte Carlo moves are wasted on them.

### Base-pair jumps
* `core/import_pose/RNA_JumpMover.cc`: a jump drawn from the D-RNA jump library between two L-RNA
  residues is mirrored. This covers the initial pairings and later `random_jump_change` moves.

### Pucker and ideal geometry
* `core/pose/rna/util.cc`: `assign_pucker()`, `apply_pucker()` (non-idealized branch), and
  `apply_ideal_c2endo_sugar_coords()` (mirrored torsions and sugar coordinates).
* `core/pose/rna/RNA_IdealCoord.{hh,cc}`: L residues (`na_lra` etc.) get mirrored reference
  geometry. Pucker detection for `ANY_PUCKER` is chirality-aware.

### Scoring (every RNA term that depends on handedness)
| Term | Change |
|---|---|
| `rna_torsion` | already chirality-aware upstream (mixed D–L bonds score 0, with a warning) |
| `rna_sugar_close` | pucker classification uses the mirrored δ |
| `rna_suite`, `suiteness_bonus` | L:L suites are scored on mirrored torsions, and derivatives are mirrored back. Mixed D:L suites are skipped |
| `rna_jr_suite`, `RNA_SuiteName` | same rules as `rna_suite` |
| `fa_stack`, `stack_elec` | use chirality-normalized base frames |
| low-res `rna_base_pair`, `rna_base_axis`, `rna_base_stagger`, `rna_base_stack`, `rna_base_stack_axis`, `rna_base_backbone`, and the DMS low-res potential | use chirality-normalized base frames (through `RNA_CentroidInfo::base_stubs()`) |
| Base-pair classification (`RNA_BasePairClassifier`) | uses `base_stubs()`, so L:L edges and cis/trans are assigned like their D mirror |
| `fa_atr`, `fa_rep`, `fa_intra_rep`, `lk_nonpolar`, `geom_sol*`, `hbond_*`, `fa_elec*`, `rna_vdw`, `rna_rg`, `ch_bond` | depend only on distances and angles, so they are mirror-invariant and need no change |
| `ref` | already maps `na_lra` etc. to the D reference energies |

The derivatives of the frame-based terms stay consistent: the normalized frame only replaces `z`
with `−z`, which is still a vector fixed to the base.

### Packing and residue types
* `core/pack/rotamer_set/rna_rotamer_building_functions.cc`: χ rotamer centers, the purine syn
  check, and the 3'-phosphate ε samples are mirrored for L residues.
* `core/chemical/PatchOperation.cc`: the `L_RNA` patch now also mirrors proton-χ (2'-OH) samples,
  and gives the L types one-letter codes `1`–`4`.
* `core/chemical/AA.cc`: one-letter codes `1`–`4` and AA names `LRA`, `LRC`, `LRG`, `LUR`.

## Limitations (still D-only or unsupported)

* **Mixed D–L covalent steps.** A D residue bonded to an L residue in the same chain gets no
  fragments across the junction, no `rna_torsion`/`rna_suite` score for that step, and no
  suite name. Put D and L strands in separate chains.
* **Heterochiral (D:L) base pairs or stacks.** No data exists for these. Low-res base-pair and stacking
  terms score them with the normalized frames, which is arbitrary. Jumps between a D and an L residue
  use the unmodified D library. Treat any D:L pairing from these terms as uninformed. Use
  constraints if you know the interface.
* **RNA–protein low-res terms** (`rnp_*`, `HackAro`) are not mirrored. They are D-RNA/L-protein
  statistics, and an L-RNA/L-protein interface has no mirror counterpart in them.
* **`RNA_HelixAssembler`** (`rna_helix`, `-dump_stems`) builds D A-form helices only. For L stems,
  supply your own L helix as a chunk (`-s`). For example, mirror a D helix by negating x in the PDB and renaming
  residues to ` 0A`/` 0C`/` 0G`/` 0U`.
* **Stepwise (SWM/SWA) samplers** and other non-FARFAR RNA samplers that hard-code torsion bins
  were not audited.
* **Chemical shift** terms (`rna_chem_shift`) use cos²-type geometry and are mirror-invariant in
  principle, but this has not been tested.

## Custom L-RNA params

The four custom files (`L_RAD_n.params` etc.) are exact mirrors of the stock `rna_phenix` D params,
so they duplicate the built-in `LRAD`/`LRCY`/`LRGU`/`LURA`. Corrected copies are kept, **unregistered**, in
`database/chemical/residue_type_sets/fa_standard/residue_types/nucleic/rna_l_custom/`. Fixes:
`AA LRA|LRC|LRG|LUR` (was the D code), `NA_ANALOGUE RCY` for C (was `RAD`), and
`IO_STRING ␣0X n` with one-letter codes `1`–`4`. Do **not** register them in `residue_types.txt`
alongside the `L_RNA` patch: two types would claim the same AA and one-letter code. To use them,
remove `patches/nucleic/rna/L_RNA.txt` from `patches.txt` first.

## FullModelParameters / FASTA notes

Earlier local changes (bracket handling in `FullModelParameters`, `name1()` remapping in
`check_full_model_info_OK`) are not needed here. Upstream already strips brackets and fills
`non_standard_residue_map` when it parses the FASTA (`core::import_pose::get_sequence_information`,
`core::sequence::parse_out_non_standard_residues`). The L types now have real one-letter codes, so
the `name1()` comparison succeeds. Bracketed names such as `3[LRGU]` also still work.

## Tests

`test/core/fragment/rna/L_RNA_Chirality.cxxtest.hh` (`fragment/rna` in `core.test.settings`):
* `1234` gives L-RNA types,
* mirrored torsions give an exact mirror image (all internal distances equal),
* puckers, C2'-endo ideal geometry and `RNA_IdealCoord` agree between D and mirrored L,
* fragment insertion negates torsions for L, leaves D untouched, and skips mixed windows,
* every hi-res (`rna_res_level_energy4` + suite terms) and low-res (`rna_lores`) score term is equal
  for a D pose and its L mirror.
