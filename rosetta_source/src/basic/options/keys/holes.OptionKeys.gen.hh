// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/holes.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_holes_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_holes_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace holes { extern BooleanOptionKey const holes; }
namespace holes { extern FileOptionKey const dalphaball; }
namespace holes { extern FileOptionKey const params; }
namespace holes { extern IntegerOptionKey const h_mode; }
namespace holes { extern BooleanOptionKey const water; }
namespace holes { extern BooleanOptionKey const make_pdb; }
namespace holes { extern BooleanOptionKey const make_voids; }
namespace holes { extern BooleanOptionKey const atom_scores; }
namespace holes { extern BooleanOptionKey const residue_scores; }
namespace holes { extern RealOptionKey const cav_shrink; }
namespace holes { extern StringOptionKey const minimize; }
namespace holes { extern BooleanOptionKey const debug; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
