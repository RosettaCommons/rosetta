// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/orbitals.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_orbitals_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_orbitals_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace orbitals { extern BooleanOptionKey const orbitals; }
namespace orbitals { extern BooleanOptionKey const Hpol; }
namespace orbitals { extern BooleanOptionKey const Haro; }
namespace orbitals { extern BooleanOptionKey const orbital; }
namespace orbitals { extern RealOptionKey const atomic_cutoff; }
namespace orbitals { extern BooleanOptionKey const verbose; }
namespace orbitals { extern BooleanOptionKey const bb_stats; }
namespace orbitals { extern BooleanOptionKey const sc_stats; }
namespace orbitals { extern BooleanOptionKey const orb_orb_stats; }
namespace orbitals { extern BooleanOptionKey const threeA; }
namespace orbitals { extern BooleanOptionKey const fiveA; }
namespace orbitals { extern BooleanOptionKey const no_pos_bonus; }
namespace orbitals { extern BooleanOptionKey const sc_bb; }
namespace orbitals { extern BooleanOptionKey const ten; }
namespace orbitals { extern BooleanOptionKey const dist; }
namespace orbitals { extern RealOptionKey const nbr_distance_squared; }
namespace orbitals { extern BooleanOptionKey const bb_smooth; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
