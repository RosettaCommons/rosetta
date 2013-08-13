// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/loops.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_loops_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_loops_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace loops { extern BooleanOptionKey const loops; }
namespace loops { extern StringOptionKey const cen_weights; }
namespace loops { extern StringOptionKey const cen_patch; }
namespace loops { extern FileOptionKey const input_pdb; }
namespace loops { extern StringVectorOptionKey const loop_file; }
namespace loops { extern FileOptionKey const extended_loop_file; }
namespace loops { extern FileOptionKey const mm_loop_file; }
namespace loops { extern BooleanOptionKey const fix_natsc; }
namespace loops { extern BooleanOptionKey const refine_only; }
namespace loops { extern BooleanOptionKey const fa_input; }
namespace loops { extern BooleanOptionKey const fast; }
namespace loops { extern FileOptionKey const vall_file; }
namespace loops { extern IntegerVectorOptionKey const frag_sizes; }
namespace loops { extern FileVectorOptionKey const frag_files; }
namespace loops { extern FileOptionKey const output_pdb; }
namespace loops { extern BooleanOptionKey const debug; }
namespace loops { extern BooleanOptionKey const build_initial; }
namespace loops { extern BooleanOptionKey const extended; }
namespace loops { extern BooleanOptionKey const remove_extended_loops; }
namespace loops { extern BooleanOptionKey const idealize_after_loop_close; }
namespace loops { extern BooleanOptionKey const idealize_before_loop_close; }
namespace loops { extern IntegerOptionKey const select_best_loop_from; }
namespace loops { extern IntegerOptionKey const build_attempts; }
namespace loops { extern IntegerOptionKey const grow_attempts; }
namespace loops { extern RealOptionKey const random_grow_loops_by; }
namespace loops { extern BooleanOptionKey const accept_aborted_loops; }
namespace loops { extern BooleanOptionKey const strict_loops; }
namespace loops { extern BooleanOptionKey const superimpose_native; }
namespace loops { extern IntegerVectorOptionKey const build_specific_loops; }
namespace loops { extern BooleanOptionKey const random_order; }
namespace loops { extern BooleanOptionKey const build_all_loops; }
namespace loops { extern BooleanOptionKey const fa_closure_protocol; }
namespace loops { extern RealOptionKey const combine_rate; }
namespace loops { extern StringOptionKey const remodel; }
namespace loops { extern StringOptionKey const intermedrelax; }
namespace loops { extern StringOptionKey const refine; }
namespace loops { extern StringOptionKey const relax; }
namespace loops { extern IntegerOptionKey const n_rebuild_tries; }
namespace loops { extern BooleanOptionKey const final_clean_fastrelax; }
namespace loops { extern BooleanOptionKey const relax_with_foldtree; }
namespace loops { extern RealOptionKey const constrain_rigid_segments; }
namespace loops { extern StringOptionKey const loopscores; }
namespace loops { extern BooleanOptionKey const timer; }
namespace loops { extern BooleanOptionKey const vicinity_sampling; }
namespace loops { extern RealOptionKey const vicinity_degree; }
namespace loops { extern RealOptionKey const neighbor_dist; }
namespace loops { extern IntegerOptionKey const kic_max_seglen; }
namespace loops { extern BooleanOptionKey const kic_recover_last; }
namespace loops { extern BooleanOptionKey const kic_min_after_repack; }
namespace loops { extern BooleanOptionKey const optimize_only_kic_region_sidechains_after_move; }
namespace loops { extern IntegerOptionKey const max_kic_build_attempts; }
namespace loops { extern IntegerOptionKey const remodel_kic_attempts; }
namespace loops { extern IntegerOptionKey const max_kic_perturber_samples; }
namespace loops { extern BooleanOptionKey const nonpivot_torsion_sampling; }
namespace loops { extern BooleanOptionKey const fix_ca_bond_angles; }
namespace loops { extern BooleanOptionKey const kic_use_linear_chainbreak; }
namespace loops { extern BooleanOptionKey const sample_omega_at_pre_prolines; }
namespace loops { extern BooleanOptionKey const allow_omega_move; }
namespace loops { extern BooleanOptionKey const kic_with_cartmin; }
namespace loops { extern BooleanOptionKey const allow_takeoff_torsion_move; }
namespace loops { extern IntegerOptionKey const extend_length; }
namespace loops { extern IntegerOptionKey const outer_cycles; }
namespace loops { extern IntegerOptionKey const max_inner_cycles; }
namespace loops { extern IntegerOptionKey const repack_period; }
namespace loops { extern BooleanOptionKey const repack_never; }
namespace loops { extern RealOptionKey const remodel_init_temp; }
namespace loops { extern RealOptionKey const remodel_final_temp; }
namespace loops { extern RealOptionKey const refine_init_temp; }
namespace loops { extern RealOptionKey const refine_final_temp; }
namespace loops { extern IntegerOptionKey const gapspan; }
namespace loops { extern IntegerOptionKey const spread; }
namespace loops { extern IntegerOptionKey const kinematic_wrapper_cycles; }
namespace loops { extern IntegerOptionKey const kic_num_rotamer_trials; }
namespace loops { extern BooleanOptionKey const kic_omega_sampling; }
namespace loops { extern RealOptionKey const kic_bump_overlap_factor; }
namespace loops { extern StringOptionKey const kic_cen_weights; }
namespace loops { extern StringOptionKey const kic_cen_patch; }
namespace loops { extern StringOptionKey const restrict_kic_sampling_to_torsion_string; }
namespace loops { extern BooleanOptionKey const derive_torsion_string_from_native_pose; }
namespace loops { extern BooleanOptionKey const always_remodel_full_loop; }
namespace loops { extern BooleanOptionKey const taboo_sampling; }
namespace loops { extern BooleanOptionKey const taboo_in_fa; }
namespace loops { extern BooleanOptionKey const ramp_fa_rep; }
namespace loops { extern BooleanOptionKey const ramp_rama; }
namespace loops { extern BooleanOptionKey const kic_rama2b; }
namespace loops { extern BooleanOptionKey const kic_no_centroid_min; }
namespace loops { extern BooleanOptionKey const kic_leave_centroid_after_initial_closure; }
namespace loops { extern BooleanOptionKey const kic_repack_neighbors_only; }
namespace loops { extern BooleanOptionKey const legacy_kic; }
namespace loops { extern BooleanOptionKey const alternative_closure_protocol; }
namespace loops { extern RealOptionKey const chainbreak_max_accept; }
namespace loops { extern BooleanOptionKey const debug_loop_closure; }
namespace loops { extern BooleanOptionKey const non_ideal_loop_closing; }
namespace loops { extern RealOptionKey const scored_frag_cycles; }
namespace loops { extern RealOptionKey const short_frag_cycles; }
namespace loops { extern RealOptionKey const rmsd_tol; }
namespace loops { extern RealOptionKey const chain_break_tol; }
namespace loops { extern BooleanOptionKey const random_loop; }
namespace loops { extern FileVectorOptionKey const stealfrags; }
namespace loops { extern IntegerOptionKey const stealfrags_times; }
namespace loops { extern RealOptionKey const coord_cst; }
namespace loops { extern RealOptionKey const skip_1mers; }
namespace loops { extern RealOptionKey const skip_3mers; }
namespace loops { extern RealOptionKey const skip_9mers; }
namespace loops { extern BooleanOptionKey const loop_model; }
namespace loops { extern RealOptionKey const score_filter_cutoff; }
namespace loops { extern BooleanOptionKey const loop_farlx; }
namespace loops { extern BooleanOptionKey const ccd_closure; }
namespace loops { extern BooleanOptionKey const skip_ccd_moves; }
namespace loops { extern BooleanOptionKey const no_randomize_loop; }
namespace loops { extern BooleanOptionKey const loops_subset; }
namespace loops { extern IntegerOptionKey const num_desired_loops; }
namespace loops { extern RealOptionKey const loop_combine_rate; }
namespace loops { extern RealOptionKey const final_score_filter; }
namespace loops { extern BooleanOptionKey const no_combine_if_fail; }
namespace loops { extern BooleanOptionKey const shorten_long_terminal_loop; }
namespace loops { extern IntegerOptionKey const backrub_trials; }
namespace loops { extern RealOptionKey const looprlx_cycle_ratio; }
namespace loops { extern RealOptionKey const extended_beta; }
namespace loops { extern BooleanOptionKey const shortrelax; }
namespace loops { extern BooleanOptionKey const fastrelax; }
namespace loops { extern BooleanOptionKey const no_looprebuild; }
namespace loops { extern BooleanOptionKey const allow_lig_move; }
namespace loops { extern FileOptionKey const keep_natro; }
namespace loops { extern IntegerOptionKey const refine_design_iterations; }
namespace loops { namespace loop_closure { extern BooleanOptionKey const loop_closure; } }
namespace loops { namespace loop_closure { extern StringOptionKey const loop_insert; } }
namespace loops { namespace loop_closure { extern StringOptionKey const loop_insert_rclrc; } }
namespace loops { namespace loop_closure { extern StringOptionKey const blueprint; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
