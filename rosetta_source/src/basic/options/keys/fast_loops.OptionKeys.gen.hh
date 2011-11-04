// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/fast_loops.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_fast_loops_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_fast_loops_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace fast_loops { extern BooleanOptionKey const fast_loops; }
namespace fast_loops { extern RealOptionKey const window_accept_ratio; }
namespace fast_loops { extern IntegerOptionKey const nr_scored_sampling_passes; }
namespace fast_loops { extern IntegerOptionKey const nr_scored_fragments; }
namespace fast_loops { extern IntegerOptionKey const min_breakout_good_loops; }
namespace fast_loops { extern IntegerOptionKey const min_breakout_fast_loops; }
namespace fast_loops { extern IntegerOptionKey const min_good_loops; }
namespace fast_loops { extern IntegerOptionKey const min_fast_loops; }
namespace fast_loops { extern RealOptionKey const vdw_delta; }
namespace fast_loops { extern IntegerOptionKey const give_up; }
namespace fast_loops { extern RealOptionKey const chainbreak_max; }
namespace fast_loops { extern FileOptionKey const fragsample_score; }
namespace fast_loops { extern FileOptionKey const fragsample_patch; }
namespace fast_loops { extern FileOptionKey const overwrite_filter_scorefxn; }
namespace fast_loops { extern FileOptionKey const patch_filter_scorefxn; }
namespace fast_loops { extern FileOptionKey const filter_cst_file; }
namespace fast_loops { extern RealOptionKey const filter_cst_weight; }
namespace fast_loops { extern FileOptionKey const fast_relax_sequence_file; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
