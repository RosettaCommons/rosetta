// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/loopfcst.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_loopfcst_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_loopfcst_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace loopfcst { extern BooleanOptionKey const loopfcst; }
namespace loopfcst { extern RealOptionKey const coord_cst_weight; }
namespace loopfcst { extern BooleanOptionKey const coord_cst_all_atom; }
namespace loopfcst { extern BooleanOptionKey const use_general_protocol; }
namespace loopfcst { extern FileOptionKey const coord_cst_weight_array; }
namespace loopfcst { extern FileOptionKey const dump_coord_cst_weight_array; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
