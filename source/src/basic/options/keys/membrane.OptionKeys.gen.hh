// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/membrane.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_membrane_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_membrane_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace membrane { extern BooleanOptionKey const membrane; }
namespace membrane { extern FileVectorOptionKey const lipid_acc_files; }
namespace membrane { extern FileVectorOptionKey const span_files; }
namespace membrane { extern FileVectorOptionKey const embed_files; }
namespace membrane { extern BooleanOptionKey const include_lips; }
namespace membrane { extern IntegerOptionKey const normal_cycles; }
namespace membrane { extern RealOptionKey const normal_mag; }
namespace membrane { extern RealOptionKey const center_mag; }
namespace membrane { extern RealOptionKey const smooth_move_frac; }
namespace membrane { extern BooleanOptionKey const no_interpolate_Mpair; }
namespace membrane { extern BooleanOptionKey const Menv_penalties; }
namespace membrane { extern BooleanOptionKey const Membed_init; }
namespace membrane { extern BooleanOptionKey const Fa_Membed_update; }
namespace membrane { extern BooleanOptionKey const center_search; }
namespace membrane { extern BooleanOptionKey const normal_search; }
namespace membrane { extern IntegerOptionKey const center_max_delta; }
namespace membrane { extern IntegerOptionKey const normal_start_angle; }
namespace membrane { extern IntegerOptionKey const normal_delta_angle; }
namespace membrane { extern IntegerOptionKey const normal_max_angle; }
namespace membrane { extern BooleanOptionKey const debug; }
namespace membrane { extern BooleanOptionKey const fixed_membrane; }
namespace membrane { extern RealVectorOptionKey const membrane_center; }
namespace membrane { extern RealVectorOptionKey const membrane_normal; }
namespace membrane { extern BooleanOptionKey const view; }
namespace membrane { extern BooleanOptionKey const Mhbond_depth; }
namespace membrane { extern RealOptionKey const thickness; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
