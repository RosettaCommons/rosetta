// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/fold_cst.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_fold_cst_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_fold_cst_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace fold_cst { extern BooleanOptionKey const fold_cst; }
namespace fold_cst { extern RealOptionKey const constraint_skip_rate; }
namespace fold_cst { extern IntegerOptionKey const violation_skip_basis; }
namespace fold_cst { extern IntegerOptionKey const violation_skip_ignore; }
namespace fold_cst { extern BooleanOptionKey const keep_skipped_csts; }
namespace fold_cst { extern BooleanOptionKey const no_minimize; }
namespace fold_cst { extern BooleanOptionKey const force_minimize; }
namespace fold_cst { extern RealVectorOptionKey const seq_sep_stages; }
namespace fold_cst { extern IntegerOptionKey const reramp_cst_cycles; }
namespace fold_cst { extern RealOptionKey const reramp_start_cstweight; }
namespace fold_cst { extern IntegerOptionKey const reramp_iterations; }
namespace fold_cst { extern BooleanOptionKey const skip_on_noviolation_in_stage1; }
namespace fold_cst { extern RealOptionKey const stage1_ramp_cst_cycle_factor; }
namespace fold_cst { extern RealOptionKey const stage2_constraint_threshold; }
namespace fold_cst { extern BooleanOptionKey const ignore_sequence_seperation; }
namespace fold_cst { extern BooleanOptionKey const no_recover_low_at_constraint_switch; }
namespace fold_cst { extern BooleanOptionKey const ramp_coord_cst; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
