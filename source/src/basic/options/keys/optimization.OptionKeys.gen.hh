// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/optimization.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_optimization_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_optimization_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace optimization { extern BooleanOptionKey const optimization; }
namespace optimization { extern IntegerOptionKey const default_max_cycles; }
namespace optimization { extern RealOptionKey const armijo_min_stepsize; }
namespace optimization { extern RealOptionKey const scale_normalmode_dampen; }
namespace optimization { extern IntegerOptionKey const lbfgs_M; }
namespace optimization { extern RealOptionKey const scale_d; }
namespace optimization { extern RealOptionKey const scale_theta; }
namespace optimization { extern RealOptionKey const scale_rb; }
namespace optimization { extern RealOptionKey const scale_rbangle; }
namespace optimization { extern BooleanOptionKey const scmin_nonideal; }
namespace optimization { extern BooleanOptionKey const scmin_cartesian; }
namespace optimization { extern BooleanOptionKey const nonideal; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
