// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/wum.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_wum_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_wum_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace wum { extern BooleanOptionKey const wum; }
namespace wum { extern IntegerOptionKey const n_slaves_per_master; }
namespace wum { extern IntegerOptionKey const n_masters; }
namespace wum { extern IntegerOptionKey const memory_limit; }
namespace wum { extern StringOptionKey const extra_scorefxn; }
namespace wum { extern FileOptionKey const extra_scorefxn_ref_structure; }
namespace wum { extern IntegerOptionKey const extra_scorefxn_relax; }
namespace wum { extern RealOptionKey const trim_proportion; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
