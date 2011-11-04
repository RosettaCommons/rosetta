// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/constraints.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_constraints_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_constraints_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace constraints { extern BooleanOptionKey const constraints; }
namespace constraints { extern RealOptionKey const CA_tether; }
namespace constraints { extern BooleanOptionKey const exit_on_bad_read; }
namespace constraints { extern StringVectorOptionKey const cst_file; }
namespace constraints { extern RealOptionKey const cst_weight; }
namespace constraints { extern RealOptionKey const max_cst_dist; }
namespace constraints { extern StringVectorOptionKey const cst_fa_file; }
namespace constraints { extern RealOptionKey const cst_fa_weight; }
namespace constraints { extern BooleanOptionKey const normalize_mixture_func; }
namespace constraints { extern BooleanOptionKey const penalize_mixture_func; }
namespace constraints { extern FileOptionKey const forest_file; }
namespace constraints { extern BooleanOptionKey const compute_total_dist_cst; }
namespace constraints { extern IntegerOptionKey const cull_with_native; }
namespace constraints { extern FileOptionKey const dump_cst_set; }
namespace constraints { extern IntegerVectorOptionKey const evaluate_max_seq_sep; }
namespace constraints { extern BooleanOptionKey const named; }
namespace constraints { extern BooleanOptionKey const no_cst_in_relax; }
namespace constraints { extern BooleanOptionKey const no_linearize_bounded; }
namespace constraints { extern RealOptionKey const pocket_constraint_weight; }
namespace constraints { extern BooleanOptionKey const pocket_zero_derivatives; }
namespace constraints { extern BooleanOptionKey const viol; }
namespace constraints { extern IntegerOptionKey const viol_level; }
namespace constraints { extern StringOptionKey const viol_type; }
namespace constraints { extern RealOptionKey const sog_cst_param; }
namespace constraints { extern BooleanOptionKey const epr_distance; }
namespace constraints { extern IntegerOptionKey const combine; }
namespace constraints { extern FileOptionKey const combine_exclude_region; }
namespace constraints { extern BooleanOptionKey const skip_redundant; }
namespace constraints { extern IntegerOptionKey const skip_redundant_width; }
namespace constraints { extern RealOptionKey const increase_constraints; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
