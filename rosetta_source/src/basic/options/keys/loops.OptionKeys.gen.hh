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

#include <platform/types.hh>
#include <utility/down_cast.hh>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
#include <utility/keys/KeyLess.fwd.hh>
#include <utility/keys/KeyLookup.fwd.hh>
#include <utility/keys/KeyLookup.hh>
#include <utility/keys/NoClient.fwd.hh>
#include <utility/keys/NoClient.hh>
#include <utility/keys/SmallKeyVector.fwd.hh>
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/AnyOption.fwd.hh>
#include <utility/options/BooleanOption.fwd.hh>
#include <utility/options/BooleanVectorOption.fwd.hh>
#include <utility/options/FileOption.fwd.hh>
#include <utility/options/FileVectorOption.fwd.hh>
#include <utility/options/IntegerOption.fwd.hh>
#include <utility/options/IntegerVectorOption.fwd.hh>
#include <utility/options/Option.fwd.hh>
#include <utility/options/PathOption.fwd.hh>
#include <utility/options/PathVectorOption.fwd.hh>
#include <utility/options/RealOption.fwd.hh>
#include <utility/options/RealVectorOption.fwd.hh>
#include <utility/options/StringOption.fwd.hh>
#include <utility/options/StringVectorOption.fwd.hh>
#include <utility/options/keys/AnyOptionKey.fwd.hh>
#include <utility/options/keys/AnyOptionKey.hh>
#include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
#include <utility/options/keys/AnyVectorOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.fwd.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
#include <utility/options/keys/BooleanVectorOptionKey.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/options/keys/IntegerOptionKey.fwd.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/OptionKeys.hh>
#include <utility/options/keys/PathOptionKey.fwd.hh>
#include <utility/options/keys/PathOptionKey.hh>
#include <utility/options/keys/PathVectorOptionKey.fwd.hh>
#include <utility/options/keys/PathVectorOptionKey.hh>
#include <utility/options/keys/RealOptionKey.fwd.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/RealVectorOptionKey.fwd.hh>
#include <utility/options/keys/RealVectorOptionKey.hh>
#include <utility/options/keys/ScalarOptionKey.fwd.hh>
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/keys/StringOptionKey.fwd.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/StringVectorOptionKey.fwd.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>
#include <utility/options/keys/VectorOptionKey.fwd.hh>
#include <utility/options/keys/VectorOptionKey.hh>
#include <utility/options/keys/all.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <utility>


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
namespace loops { extern BooleanOptionKey const allow_omega_move; }
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

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
