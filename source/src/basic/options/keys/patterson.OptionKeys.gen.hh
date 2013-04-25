// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/patterson.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_patterson_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_patterson_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace patterson { extern BooleanOptionKey const patterson; }
namespace patterson { extern BooleanOptionKey const debug; }
namespace patterson { extern RealOptionKey const weight; }
namespace patterson { extern RealOptionKey const sc_scaling; }
namespace patterson { extern RealVectorOptionKey const radius_cutoffs; }
namespace patterson { extern RealVectorOptionKey const resolution_cutoffs; }
namespace patterson { extern RealOptionKey const Bsol; }
namespace patterson { extern RealOptionKey const Fsol; }
namespace patterson { extern RealOptionKey const model_B; }
namespace patterson { extern RealOptionKey const rmsd; }
namespace patterson { extern BooleanOptionKey const no_ecalc; }
namespace patterson { extern IntegerOptionKey const nshells; }
namespace patterson { extern BooleanOptionKey const use_spline_interpolation; }
namespace patterson { extern BooleanOptionKey const use_on_repack; }
namespace patterson { extern BooleanOptionKey const dont_use_symm_in_pcalc; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
