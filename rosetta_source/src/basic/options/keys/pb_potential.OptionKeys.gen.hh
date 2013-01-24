// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/pb_potential.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_pb_potential_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_pb_potential_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace pb_potential { extern BooleanOptionKey const pb_potential; }
namespace pb_potential { extern IntegerVectorOptionKey const charged_chains; }
namespace pb_potential { extern BooleanOptionKey const sidechain_only; }
namespace pb_potential { extern IntegerVectorOptionKey const revamp_near_chain; }
namespace pb_potential { extern StringOptionKey const apbs_path; }
namespace pb_potential { extern RealOptionKey const potential_cap; }
namespace pb_potential { extern RealOptionKey const epsilon; }
namespace pb_potential { extern IntegerOptionKey const apbs_debug; }
namespace pb_potential { extern BooleanOptionKey const calcenergy; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
