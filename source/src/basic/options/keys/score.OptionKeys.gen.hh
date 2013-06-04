// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/score.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_score_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_score_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace score { extern BooleanOptionKey const score_pose_cutpoint_variants; }
namespace score { extern BooleanOptionKey const score; }
namespace score { extern StringOptionKey const weights; }
namespace score { extern StringVectorOptionKey const set_weights; }
namespace score { extern StringOptionKey const pack_weights; }
namespace score { extern StringOptionKey const soft_wts; }
namespace score { extern BooleanOptionKey const docking_interface_score; }
namespace score { extern RealOptionKey const min_score_score; }
namespace score { extern StringOptionKey const custom_atom_pair; }
namespace score { extern FileVectorOptionKey const patch; }
namespace score { extern BooleanOptionKey const empty; }
namespace score { extern RealOptionKey const fa_max_dis; }
namespace score { extern BooleanOptionKey const fa_Hatr; }
namespace score { extern BooleanOptionKey const no_smooth_etables; }
namespace score { extern RealOptionKey const etable_lr; }
namespace score { extern BooleanOptionKey const no_lk_polar_desolvation; }
namespace score { extern StringOptionKey const input_etables; }
namespace score { extern StringOptionKey const output_etables; }
namespace score { extern BooleanOptionKey const analytic_etable_evaluation; }
namespace score { extern RealOptionKey const rms_target; }
namespace score { extern BooleanOptionKey const ramaneighbors; }
namespace score { extern StringOptionKey const optH_weights; }
namespace score { extern StringOptionKey const optH_patch; }
namespace score { extern StringOptionKey const hbond_params; }
namespace score { extern BooleanOptionKey const hbond_disable_bbsc_exclusion_rule; }
namespace score { extern IntegerOptionKey const symE_units; }
namespace score { extern RealOptionKey const symE_bonus; }
namespace score { extern RealOptionKey const NV_lbound; }
namespace score { extern RealOptionKey const NV_ubound; }
namespace score { extern StringOptionKey const NV_table; }
namespace score { extern BooleanOptionKey const disable_orientation_dependent_rna_ch_o_bonds; }
namespace score { extern StringOptionKey const rna_torsion_potential; }
namespace score { extern BooleanOptionKey const rna_torsion_skip_chainbreak; }
namespace score { extern StringOptionKey const rna_chemical_shift_exp_data; }
namespace score { extern StringOptionKey const rna_chemical_shift_H5_prime_mode; }
namespace score { extern IntegerVectorOptionKey const rna_chemical_shift_include_res; }
namespace score { extern BooleanOptionKey const use_2prime_OH_potential; }
namespace score { extern BooleanOptionKey const include_neighbor_base_stacks; }
namespace score { extern BooleanOptionKey const find_neighbors_3dgrid; }
namespace score { extern BooleanOptionKey const find_neighbors_stripehash; }
namespace score { extern StringOptionKey const seqdep_refene_fname; }
namespace score { extern StringOptionKey const secondary_seqdep_refene_fname; }
namespace score { extern BooleanOptionKey const exact_occ_pairwise; }
namespace score { extern BooleanOptionKey const exact_occ_skip_Hbonders; }
namespace score { extern BooleanOptionKey const exact_occ_include_Hbond_contribution; }
namespace score { extern BooleanOptionKey const exact_occ_pairwise_by_res; }
namespace score { extern BooleanOptionKey const exact_occ_split_between_res; }
namespace score { extern BooleanOptionKey const exact_occ_self_res_no_occ; }
namespace score { extern RealOptionKey const exact_occ_radius_scaling; }
namespace score { extern StringVectorOptionKey const ref_offsets; }
namespace score { extern BooleanOptionKey const output_residue_energies; }
namespace score { extern StringOptionKey const fa_custom_pair_distance_file; }
namespace score { extern RealOptionKey const disulf_matching_probe; }
namespace score { extern RealVectorOptionKey const bonded_params; }
namespace score { extern StringOptionKey const bonded_params_dir; }
namespace score { extern StringOptionKey const extra_improper_file; }
namespace score { extern RealOptionKey const pro_close_planar_constraint; }
namespace score { extern BooleanOptionKey const linear_bonded_potential; }
namespace score { extern BooleanOptionKey const geom_sol_correct_acceptor_base; }
namespace score { extern IntegerVectorOptionKey const rg_local_span; }
namespace score { extern BooleanOptionKey const unmodifypot; }
namespace score { namespace saxs { extern BooleanOptionKey const saxs; } }
namespace score { namespace saxs { extern RealOptionKey const min_score; } }
namespace score { namespace saxs { extern StringOptionKey const custom_ff; } }
namespace score { namespace saxs { extern StringOptionKey const print_i_calc; } }
namespace score { namespace saxs { extern FileOptionKey const ref_fa_spectrum; } }
namespace score { namespace saxs { extern FileOptionKey const ref_cen_spectrum; } }
namespace score { namespace saxs { extern FileOptionKey const ref_spectrum; } }
namespace score { namespace saxs { extern FileOptionKey const ref_pddf; } }
namespace score { namespace saxs { extern BooleanOptionKey const skip_hydrogens; } }
namespace score { namespace saxs { extern RealOptionKey const d_min; } }
namespace score { namespace saxs { extern RealOptionKey const d_max; } }
namespace score { namespace saxs { extern RealOptionKey const d_step; } }
namespace score { namespace saxs { extern RealOptionKey const q_min; } }
namespace score { namespace saxs { extern RealOptionKey const q_max; } }
namespace score { namespace saxs { extern RealOptionKey const q_step; } }
namespace score { namespace saxs { extern BooleanOptionKey const fit_pddf_area; } }
namespace score { extern IntegerVectorOptionKey const sidechain_buried; }
namespace score { extern IntegerVectorOptionKey const sidechain_exposed; }
namespace score { extern RealOptionKey const hackelec_min_dis; }
namespace score { extern RealOptionKey const hackelec_max_dis; }
namespace score { extern RealOptionKey const hackelec_die; }
namespace score { extern BooleanOptionKey const hackelec_r_option; }
namespace score { extern BooleanOptionKey const smooth_hack_elec; }
namespace score { extern RealOptionKey const facts_GBpair_cut; }
namespace score { extern RealOptionKey const facts_min_dis; }
namespace score { extern RealOptionKey const facts_kappa; }
namespace score { extern BooleanOptionKey const facts_apprx; }
namespace score { extern IntegerOptionKey const facts_asp_patch; }
namespace score { extern RealOptionKey const facts_selfenergy_scale; }
namespace score { extern BooleanOptionKey const facts_plane_to_self; }
namespace score { extern RealOptionKey const facts_intrares_scale; }
namespace score { extern RealOptionKey const facts_saltbridge_correction; }
namespace score { extern RealOptionKey const facts_dshift; }
namespace score { extern RealOptionKey const facts_die; }
namespace score { extern RealOptionKey const facts_elec_sh_exponent; }
namespace score { extern StringOptionKey const nmer_ref_energies; }
namespace score { extern StringOptionKey const nmer_ref_energies_list; }
namespace score { extern StringOptionKey const nmer_pssm; }
namespace score { extern StringOptionKey const nmer_pssm_list; }
namespace score { extern RealOptionKey const nmer_pssm_scorecut; }
namespace score { extern StringOptionKey const nmer_svm; }
namespace score { extern StringOptionKey const nmer_svm_list; }
namespace score { extern RealOptionKey const nmer_svm_scorecut; }
namespace score { extern StringOptionKey const nmer_svm_aa_matrix; }
namespace score { extern IntegerOptionKey const nmer_svm_term_length; }
namespace score { extern BooleanOptionKey const nmer_svm_pssm_feat; }
namespace score { extern IntegerOptionKey const nmer_ref_seq_length; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
