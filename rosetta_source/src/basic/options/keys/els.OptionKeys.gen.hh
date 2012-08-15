// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/els.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_els_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_els_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace els { extern BooleanOptionKey const els; }
namespace els { extern IntegerOptionKey const master_wu_per_send; }
namespace els { extern StringOptionKey const vars; }
namespace els { extern FileOptionKey const script; }
namespace els { extern IntegerOptionKey const num_traj; }
namespace els { extern IntegerOptionKey const traj_per_master; }
namespace els { extern IntegerOptionKey const shortest_wu; }
namespace els { extern BooleanOptionKey const pool; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
