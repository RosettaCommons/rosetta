// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/flxbb.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_flxbb_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_flxbb_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace flxbb { extern BooleanOptionKey const flxbb; }
namespace flxbb { extern BooleanOptionKey const view; }
namespace flxbb { extern IntegerOptionKey const ncycle; }
namespace flxbb { extern RealOptionKey const constraints_sheet; }
namespace flxbb { extern BooleanOptionKey const constraints_sheet_include_cacb_pseudotorsion; }
namespace flxbb { extern RealOptionKey const constraints_NtoC; }
namespace flxbb { extern IntegerOptionKey const filter_trial; }
namespace flxbb { extern StringOptionKey const filter_type; }
namespace flxbb { extern BooleanOptionKey const exclude_Met; }
namespace flxbb { extern BooleanOptionKey const exclude_Ala; }
namespace flxbb { extern FileOptionKey const blueprint; }
namespace flxbb { extern BooleanOptionKey const movemap_from_blueprint; }
namespace flxbb { namespace layer { extern StringOptionKey const layer; } }
namespace flxbb { namespace layer { extern RealOptionKey const pore_radius; } }
namespace flxbb { namespace layer { extern RealOptionKey const burial; } }
namespace flxbb { namespace layer { extern RealOptionKey const surface; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
