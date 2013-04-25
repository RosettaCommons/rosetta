// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/run.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_run_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_run_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace run { extern BooleanOptionKey const run; }
namespace run { extern FileVectorOptionKey const batches; }
namespace run { extern BooleanOptionKey const no_prof_info_in_silentout; }
namespace run { extern BooleanOptionKey const archive; }
namespace run { extern IntegerOptionKey const n_replica; }
namespace run { extern BooleanOptionKey const shuffle; }
namespace run { extern IntegerOptionKey const n_cycles; }
namespace run { extern IntegerOptionKey const repeat; }
namespace run { extern IntegerOptionKey const max_min_iter; }
namespace run { extern IntegerOptionKey const maxruntime; }
namespace run { extern BooleanOptionKey const write_failures; }
namespace run { extern BooleanOptionKey const clean; }
namespace run { extern BooleanOptionKey const benchmark; }
namespace run { extern BooleanOptionKey const test_cycles; }
namespace run { extern BooleanOptionKey const memory_test_cycles; }
namespace run { extern BooleanOptionKey const dry_run; }
namespace run { extern BooleanOptionKey const debug; }
namespace run { extern BooleanOptionKey const profile; }
namespace run { extern IntegerOptionKey const max_retry_job; }
namespace run { extern IntegerOptionKey const verbosity; }
namespace run { extern BooleanOptionKey const version; }
namespace run { extern BooleanOptionKey const nodelay; }
namespace run { extern IntegerOptionKey const delay; }
namespace run { extern IntegerOptionKey const random_delay; }
namespace run { extern BooleanOptionKey const timer; }
namespace run { extern StringOptionKey const series; }
namespace run { extern StringOptionKey const protein; }
namespace run { extern StringOptionKey const chain; }
namespace run { extern BooleanOptionKey const score_only; }
namespace run { extern BooleanOptionKey const silent_input; }
namespace run { extern BooleanOptionKey const decoystats; }
namespace run { extern BooleanOptionKey const output_hbond_info; }
namespace run { extern RealOptionKey const wide_nblist_extension; }
namespace run { extern BooleanOptionKey const status; }
namespace run { extern BooleanOptionKey const constant_seed; }
namespace run { extern IntegerOptionKey const jran; }
namespace run { extern BooleanOptionKey const use_time_as_seed; }
namespace run { extern StringOptionKey const rng_seed_device; }
namespace run { extern IntegerOptionKey const seed_offset; }
namespace run { extern StringOptionKey const rng; }
namespace run { extern IntegerOptionKey const run_level; }
namespace run { extern StringOptionKey const verbose; }
namespace run { extern BooleanOptionKey const silent; }
namespace run { extern BooleanOptionKey const regions; }
namespace run { extern BooleanOptionKey const find_disulf; }
namespace run { extern BooleanOptionKey const rebuild_disulf; }
namespace run { extern BooleanOptionKey const movie; }
namespace run { extern BooleanOptionKey const trajectory; }
namespace run { extern BooleanOptionKey const IUPAC; }
namespace run { extern BooleanOptionKey const preserve_header; }
namespace run { extern BooleanOptionKey const evolution; }
namespace run { extern BooleanOptionKey const suppress_checkpoints; }
namespace run { extern BooleanOptionKey const checkpoint; }
namespace run { extern BooleanOptionKey const delete_checkpoints; }
namespace run { extern IntegerOptionKey const checkpoint_interval; }
namespace run { extern StringOptionKey const protocol; }
namespace run { extern BooleanOptionKey const remove_ss_length_screen; }
namespace run { extern StringOptionKey const min_type; }
namespace run { extern RealOptionKey const min_tolerance; }
namespace run { extern BooleanOptionKey const nblist_autoupdate; }
namespace run { extern RealOptionKey const nblist_autoupdate_narrow; }
namespace run { extern RealOptionKey const nblist_autoupdate_wide; }
namespace run { extern BooleanOptionKey const skip_set_reasonable_fold_tree; }
namespace run { extern BooleanOptionKey const randomize_missing_coords; }
namespace run { extern BooleanOptionKey const ignore_zero_occupancy; }
namespace run { extern IntegerOptionKey const cycles_outer; }
namespace run { extern IntegerOptionKey const cycles_inner; }
namespace run { extern IntegerOptionKey const repack_rate; }
namespace run { extern BooleanOptionKey const reinitialize_mover_for_each_job; }
namespace run { extern BooleanOptionKey const reinitialize_mover_for_new_input; }
namespace run { extern BooleanOptionKey const multiple_processes_writing_to_one_directory; }
namespace run { extern StringOptionKey const jobdist_miscfile_ext; }
namespace run { extern BooleanOptionKey const no_scorefile; }
namespace run { extern BooleanOptionKey const other_pose_to_scorefile; }
namespace run { extern FileOptionKey const other_pose_scorefile; }
namespace run { extern BooleanOptionKey const intermediate_scorefiles; }
namespace run { extern BooleanOptionKey const intermediate_structures; }
namespace run { extern BooleanOptionKey const idealize_before_protocol; }
namespace run { extern BooleanOptionKey const interactive; }
namespace run { extern BooleanOptionKey const condor; }
namespace run { extern IntegerOptionKey const nproc; }
namespace run { extern IntegerOptionKey const proc_id; }
namespace run { extern BooleanOptionKey const exit_if_missing_heavy_atoms; }
namespace run { extern RealOptionKey const show_simulation_in_pymol; }
namespace run { extern BooleanOptionKey const keep_pymol_simulation_history; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
