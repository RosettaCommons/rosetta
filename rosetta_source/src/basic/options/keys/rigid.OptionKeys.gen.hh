// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/rigid.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_rigid_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_rigid_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace rigid { extern BooleanOptionKey const rigid; }
namespace rigid { extern IntegerOptionKey const fragment_cycles; }
namespace rigid { extern IntegerOptionKey const jump_attempts; }
namespace rigid { extern IntegerOptionKey const jump_cycles; }
namespace rigid { extern RealOptionKey const max_ca_ca_dist; }
namespace rigid { extern IntegerOptionKey const rigid_body_cycles; }
namespace rigid { extern RealOptionKey const rotation; }
namespace rigid { extern IntegerOptionKey const sequence_separation; }
namespace rigid { extern IntegerOptionKey const small_cycles; }
namespace rigid { extern RealOptionKey const temperature; }
namespace rigid { extern RealOptionKey const translation; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
