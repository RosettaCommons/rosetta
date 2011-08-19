// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/abinitio.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_abinitio_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_abinitio_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace abinitio { extern RealOptionKey const prob_perturb_weights; }
namespace abinitio { extern BooleanOptionKey const abinitio; }
namespace abinitio { extern BooleanOptionKey const membrane; }
namespace abinitio { extern FileOptionKey const kill_hairpins; }
namespace abinitio { extern RealOptionKey const kill_hairpins_frequency; }
namespace abinitio { extern BooleanOptionKey const smooth_cycles_only; }
namespace abinitio { extern BooleanOptionKey const relax; }
namespace abinitio { extern BooleanOptionKey const final_clean_relax; }
namespace abinitio { extern BooleanOptionKey const fastrelax; }
namespace abinitio { extern BooleanOptionKey const multifastrelax; }
namespace abinitio { extern BooleanOptionKey const debug; }
namespace abinitio { extern BooleanOptionKey const clear_pose_cache; }
namespace abinitio { extern BooleanOptionKey const explicit_pdb_debug; }
namespace abinitio { extern BooleanOptionKey const use_filters; }
namespace abinitio { extern RealOptionKey const increase_cycles; }
namespace abinitio { extern IntegerOptionKey const number_3mer_frags; }
namespace abinitio { extern IntegerOptionKey const number_9mer_frags; }
namespace abinitio { extern RealOptionKey const temperature; }
namespace abinitio { extern RealOptionKey const rg_reweight; }
namespace abinitio { extern RealOptionKey const strand_dist_cutoff; }
namespace abinitio { extern BooleanOptionKey const stretch_strand_dist_cutoff; }
namespace abinitio { extern RealOptionKey const rsd_wt_helix; }
namespace abinitio { extern RealOptionKey const rsd_wt_strand; }
namespace abinitio { extern RealOptionKey const rsd_wt_loop; }
namespace abinitio { extern BooleanOptionKey const fast; }
namespace abinitio { extern BooleanOptionKey const skip_convergence_check; }
namespace abinitio { extern FileVectorOptionKey const stage1_patch; }
namespace abinitio { extern FileVectorOptionKey const stage2_patch; }
namespace abinitio { extern FileVectorOptionKey const stage3a_patch; }
namespace abinitio { extern FileVectorOptionKey const stage3b_patch; }
namespace abinitio { extern FileVectorOptionKey const stage4_patch; }
namespace abinitio { extern FileVectorOptionKey const stage5_patch; }
namespace abinitio { extern BooleanOptionKey const exit_when_converged; }
namespace abinitio { extern BooleanOptionKey const steal_3mers; }
namespace abinitio { extern BooleanOptionKey const steal_9mers; }
namespace abinitio { extern BooleanOptionKey const no_write_failures; }
namespace abinitio { extern BooleanOptionKey const relax_failures; }
namespace abinitio { extern BooleanOptionKey const relax_with_jumps; }
namespace abinitio { extern BooleanOptionKey const process_store; }
namespace abinitio { extern IntegerVectorOptionKey const fix_residues_to_native; }
namespace abinitio { extern BooleanOptionKey const return_full_atom; }
namespace abinitio { extern BooleanOptionKey const detect_disulfide_before_relax; }
namespace abinitio { extern BooleanOptionKey const close_loops; }
namespace abinitio { extern BooleanOptionKey const bGDT; }
namespace abinitio { extern BooleanOptionKey const dump_frags; }
namespace abinitio { extern BooleanOptionKey const jdist_rerun; }
namespace abinitio { extern RealOptionKey const perturb; }
namespace abinitio { extern BooleanOptionKey const rerun; }
namespace abinitio { extern IntegerVectorOptionKey const rmsd_residues; }
namespace abinitio { extern BooleanOptionKey const start_native; }
namespace abinitio { extern BooleanOptionKey const debug_structures; }
namespace abinitio { extern FileOptionKey const log_frags; }
namespace abinitio { extern BooleanOptionKey const only_stage1; }
namespace abinitio { extern RealOptionKey const end_bias; }
namespace abinitio { extern IntegerOptionKey const symmetry_residue; }
namespace abinitio { extern RealOptionKey const vdw_weight_stage1; }
namespace abinitio { extern BooleanOptionKey const override_vdw_all_stages; }
namespace abinitio { extern IntegerVectorOptionKey const recover_low_in_stages; }
namespace abinitio { extern IntegerVectorOptionKey const skip_stages; }
namespace abinitio { extern BooleanOptionKey const close_chbrk; }
namespace abinitio { extern BooleanOptionKey const include_stage5; }
namespace abinitio { extern BooleanOptionKey const close_loops_by_idealizing; }
namespace abinitio { extern BooleanOptionKey const optimize_cutpoints_using_kic; }
namespace abinitio { extern IntegerOptionKey const optimize_cutpoints_margin; }
namespace abinitio { extern FileOptionKey const HD_EX_Info; }
namespace abinitio { extern RealOptionKey const HD_penalty; }
namespace abinitio { extern RealOptionKey const HD_fa_penalty; }
namespace abinitio { extern FileOptionKey const sheet_edge_pred; }
namespace abinitio { extern RealOptionKey const SEP_score_scalling; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
