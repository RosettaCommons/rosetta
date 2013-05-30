// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/mistakes.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_mistakes_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_mistakes_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace mistakes { extern BooleanOptionKey const mistakes; }
namespace mistakes { extern BooleanOptionKey const restore_pre_talaris_2013_behavior; }
namespace mistakes { namespace chemical { extern BooleanOptionKey const chemical; } }
namespace mistakes { namespace chemical { extern BooleanOptionKey const pre_talaris2013_geometries; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
