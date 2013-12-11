// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/in.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_in_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_in_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace in { extern BooleanOptionKey const in; }
namespace in { extern StringOptionKey const Ntermini; }
namespace in { extern StringOptionKey const Ctermini; }
namespace in { extern BooleanOptionKey const use_truncated_termini; }
namespace in { extern BooleanOptionKey const ignore_unrecognized_res; }
namespace in { extern BooleanOptionKey const ignore_waters; }
namespace in { extern BooleanOptionKey const add_orbitals; }
namespace in { extern BooleanOptionKey const show_all_fixes; }
namespace in { extern BooleanOptionKey const include_sugars; }
namespace in { extern BooleanOptionKey const include_surfaces; }
namespace in { extern BooleanOptionKey const enable_branching; }
namespace in { extern BooleanOptionKey const remember_unrecognized_res; }
namespace in { extern BooleanOptionKey const remember_unrecognized_water; }
namespace in { extern BooleanOptionKey const preserve_crystinfo; }
namespace in { extern BooleanOptionKey const detect_oops; }
namespace in { extern BooleanOptionKey const detect_disulf; }
namespace in { extern RealOptionKey const detect_disulf_tolerance; }
namespace in { extern FileOptionKey const fix_disulf; }
namespace in { extern BooleanOptionKey const missing_density_to_jump; }
namespace in { extern BooleanOptionKey const use_stupid_foldtree_format; }
namespace in { extern IntegerVectorOptionKey const target_residues; }
namespace in { extern IntegerVectorOptionKey const replonly_residues; }
namespace in { extern BooleanOptionKey const replonly_loops; }
namespace in { extern BooleanOptionKey const use_database; }
namespace in { namespace dbms { extern BooleanOptionKey const dbms; } }
namespace in { namespace dbms { extern StringVectorOptionKey const struct_ids; } }
namespace in { extern IntegerOptionKey const database_protocol; }
namespace in { extern StringVectorOptionKey const select_structures_from_database; }
namespace in { namespace path { extern PathVectorOptionKey const path; } }
namespace in { namespace path { extern PathVectorOptionKey const fragments; } }
namespace in { namespace path { extern PathVectorOptionKey const pdb; } }
namespace in { namespace path { extern PathVectorOptionKey const database; } }
namespace in { namespace file { extern BooleanOptionKey const file; } }
namespace in { namespace file { extern FileVectorOptionKey const s; } }
namespace in { namespace file { extern FileVectorOptionKey const l; } }
namespace in { namespace file { extern FileVectorOptionKey const list; } }
namespace in { namespace file { extern FileVectorOptionKey const screening_list; } }
namespace in { namespace file { extern FileOptionKey const screening_job_file; } }
namespace in { namespace file { extern BooleanOptionKey const shuffle_screening_jobs; } }
namespace in { namespace file { extern FileOptionKey const native; } }
namespace in { namespace file { extern FileOptionKey const torsion_bin_probs; } }
namespace in { namespace file { extern FileOptionKey const PCS_frag_cst; } }
namespace in { namespace file { extern FileOptionKey const talos_phi_psi; } }
namespace in { namespace file { extern FileOptionKey const talos_cs; } }
namespace in { namespace file { extern FileOptionKey const ambig_talos_cs_A; } }
namespace in { namespace file { extern FileOptionKey const ambig_talos_cs_B; } }
namespace in { namespace file { extern IntegerVectorOptionKey const native_exclude_res; } }
namespace in { namespace file { extern StringVectorOptionKey const tags; } }
namespace in { namespace file { extern StringVectorOptionKey const user_tags; } }
namespace in { namespace file { extern FileOptionKey const tagfile; } }
namespace in { namespace file { extern FileVectorOptionKey const frag_files; } }
namespace in { namespace file { extern IntegerVectorOptionKey const frag_sizes; } }
namespace in { namespace file { extern FileVectorOptionKey const extra_res; } }
namespace in { namespace file { extern FileVectorOptionKey const extra_res_fa; } }
namespace in { namespace file { extern FileVectorOptionKey const extra_res_mol; } }
namespace in { namespace file { extern StringOptionKey const extra_res_database; } }
namespace in { namespace file { extern StringOptionKey const extra_res_pq_schema; } }
namespace in { namespace file { extern StringOptionKey const extra_res_database_mode; } }
namespace in { namespace file { extern FileOptionKey const extra_res_database_resname_list; } }
namespace in { namespace file { extern FileVectorOptionKey const extra_res_cen; } }
namespace in { namespace file { extern PathVectorOptionKey const extra_res_path; } }
namespace in { namespace file { extern PathVectorOptionKey const extra_res_batch_path; } }
namespace in { namespace file { extern FileVectorOptionKey const extra_patch_fa; } }
namespace in { namespace file { extern FileVectorOptionKey const extra_patch_cen; } }
namespace in { namespace file { extern StringOptionKey const frag3; } }
namespace in { namespace file { extern StringOptionKey const frag9; } }
namespace in { namespace file { extern StringOptionKey const fragA; } }
namespace in { namespace file { extern StringOptionKey const fragB; } }
namespace in { namespace file { extern StringOptionKey const surface_vectors; } }
namespace in { namespace file { extern StringOptionKey const xyz; } }
namespace in { namespace file { extern IntegerOptionKey const fragA_size; } }
namespace in { namespace file { extern BooleanOptionKey const keep_input_scores; } }
namespace in { namespace file { extern BooleanOptionKey const lazy_silent; } }
namespace in { namespace file { extern FileVectorOptionKey const silent; } }
namespace in { namespace file { extern FileVectorOptionKey const atom_tree_diff; } }
namespace in { namespace file { extern StringOptionKey const zip; } }
namespace in { namespace file { extern FileVectorOptionKey const boinc_wu_zip; } }
namespace in { namespace file { extern BooleanOptionKey const fullatom; } }
namespace in { namespace file { extern BooleanOptionKey const centroid_input; } }
namespace in { namespace file { extern BooleanOptionKey const centroid; } }
namespace in { namespace file { extern StringOptionKey const treat_residues_in_these_chains_as_separate_chemical_entities; } }
namespace in { namespace file { extern StringOptionKey const residue_type_set; } }
namespace in { namespace file { extern FileOptionKey const pca; } }
namespace in { namespace file { extern RealOptionKey const silent_energy_cut; } }
namespace in { namespace file { extern FileVectorOptionKey const silent_list; } }
namespace in { namespace file { extern BooleanOptionKey const silent_renumber; } }
namespace in { namespace file { extern BooleanOptionKey const silent_optH; } }
namespace in { namespace file { extern StringOptionKey const silent_struct_type; } }
namespace in { namespace file { extern BooleanOptionKey const silent_read_through_errors; } }
namespace in { namespace file { extern StringOptionKey const silent_score_prefix; } }
namespace in { namespace file { extern IntegerOptionKey const silent_select_random; } }
namespace in { namespace file { extern IntegerOptionKey const silent_select_range_start; } }
namespace in { namespace file { extern IntegerOptionKey const silent_select_range_mul; } }
namespace in { namespace file { extern IntegerOptionKey const silent_select_range_len; } }
namespace in { namespace file { extern BooleanOptionKey const skip_failed_simulations; } }
namespace in { namespace file { extern StringVectorOptionKey const silent_scores_wanted; } }
namespace in { namespace file { extern FileVectorOptionKey const fasta; } }
namespace in { namespace file { extern FileVectorOptionKey const pssm; } }
namespace in { namespace file { extern StringVectorOptionKey const seq; } }
namespace in { namespace file { extern FileOptionKey const checkpoint; } }
namespace in { namespace file { extern FileVectorOptionKey const alignment; } }
namespace in { namespace file { extern FileVectorOptionKey const alignment2; } }
namespace in { namespace file { extern FileOptionKey const rama2b_map; } }
namespace in { namespace file { extern FileOptionKey const psipred_ss2; } }
namespace in { namespace file { extern FileOptionKey const dssp; } }
namespace in { namespace file { extern BooleanOptionKey const fail_on_bad_hbond; } }
namespace in { namespace file { extern FileOptionKey const movemap; } }
namespace in { namespace file { extern BooleanOptionKey const repair_sidechains; } }
namespace in { namespace file { extern BooleanOptionKey const no_binary_dunlib; } }
namespace in { namespace file { extern IntegerOptionKey const extended_pose; } }
namespace in { namespace file { extern FileVectorOptionKey const template_pdb; } }
namespace in { namespace file { extern FileOptionKey const template_silent; } }
namespace in { namespace file { extern FileVectorOptionKey const rdc; } }
namespace in { namespace file { extern FileVectorOptionKey const csa; } }
namespace in { namespace file { extern FileVectorOptionKey const dc; } }
namespace in { namespace file { extern FileVectorOptionKey const burial; } }
namespace in { namespace file { extern FileVectorOptionKey const vall; } }
namespace in { namespace file { extern BooleanOptionKey const rescore; } }
namespace in { namespace file { extern StringOptionKey const spanfile; } }
namespace in { namespace file { extern StringOptionKey const lipofile; } }
namespace in { namespace file { extern StringOptionKey const HDX; } }
namespace in { namespace file { extern RealOptionKey const d2h_sa_reweight; } }
namespace in { namespace file { extern FileOptionKey const sucker_params; } }
namespace in { namespace file { extern FileOptionKey const fold_tree; } }
namespace in { namespace file { extern BooleanOptionKey const obey_ENDMDL; } }
namespace in { namespace file { extern BooleanOptionKey const new_chain_order; } }
namespace in { namespace file { extern FileOptionKey const ddg_predictions_file; } }
namespace in { namespace file { extern IntegerVectorOptionKey const input_res; } }
namespace in { namespace file { extern IntegerVectorOptionKey const minimize_res; } }
namespace in { namespace file { extern StringOptionKey const md_schfile; } }
namespace in { namespace file { extern BooleanOptionKey const read_pdb_link_records; } }
namespace in { namespace file { extern FileOptionKey const native_contacts; } }
namespace in { namespace rdf { extern BooleanOptionKey const rdf; } }
namespace in { namespace rdf { extern BooleanOptionKey const sep_bb_ss; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
