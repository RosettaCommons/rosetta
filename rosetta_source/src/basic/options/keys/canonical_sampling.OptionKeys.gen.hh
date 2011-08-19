// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/canonical_sampling.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_canonical_sampling_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_canonical_sampling_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace canonical_sampling { extern BooleanOptionKey const canonical_sampling; }
namespace canonical_sampling { namespace probabilities { extern BooleanOptionKey const probabilities; } }
namespace canonical_sampling { namespace probabilities { extern RealOptionKey const sc; } }
namespace canonical_sampling { namespace probabilities { extern RealOptionKey const localbb; } }
namespace canonical_sampling { namespace probabilities { extern RealOptionKey const sc_prob_uniform; } }
namespace canonical_sampling { namespace probabilities { extern RealOptionKey const sc_prob_withinrot; } }
namespace canonical_sampling { namespace probabilities { extern RealOptionKey const sc_prob_perturbcurrent; } }
namespace canonical_sampling { namespace probabilities { extern BooleanOptionKey const MPI_sync_pools; } }
namespace canonical_sampling { namespace probabilities { extern BooleanOptionKey const MPI_bcast; } }
namespace canonical_sampling { namespace probabilities { extern BooleanOptionKey const fast_sc_moves; } }
namespace canonical_sampling { namespace probabilities { extern RealOptionKey const fast_sc_moves_ntrials; } }
namespace canonical_sampling { namespace probabilities { extern BooleanOptionKey const no_jd2_output; } }
namespace canonical_sampling { namespace probabilities { extern BooleanOptionKey const use_hierarchical_clustering; } }
namespace canonical_sampling { namespace probabilities { extern IntegerOptionKey const hierarchical_max_cache_size; } }
namespace canonical_sampling { namespace probabilities { extern RealOptionKey const backrub; } }
namespace canonical_sampling { namespace probabilities { extern RealOptionKey const conrot; } }
namespace canonical_sampling { namespace sampling { extern BooleanOptionKey const sampling; } }
namespace canonical_sampling { namespace sampling { extern BooleanOptionKey const no_detailed_balance; } }
namespace canonical_sampling { namespace sampling { extern IntegerOptionKey const ntrials; } }
namespace canonical_sampling { namespace sampling { extern RealOptionKey const mc_kt; } }
namespace canonical_sampling { namespace sampling { extern IntegerOptionKey const interval_pose_dump; } }
namespace canonical_sampling { namespace sampling { extern IntegerOptionKey const interval_data_dump; } }
namespace canonical_sampling { namespace sampling { extern BooleanOptionKey const output_only_cluster_transitions; } }
namespace canonical_sampling { namespace sampling { extern RealOptionKey const transition_threshold; } }
namespace canonical_sampling { namespace sampling { extern IntegerOptionKey const max_files_per_dir; } }
namespace canonical_sampling { namespace sampling { extern BooleanOptionKey const save_loops_only; } }
namespace canonical_sampling { namespace sampling { extern BooleanOptionKey const dump_loops_only; } }
namespace canonical_sampling { namespace out { extern BooleanOptionKey const out; } }
namespace canonical_sampling { namespace out { extern FileOptionKey const new_structures; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
