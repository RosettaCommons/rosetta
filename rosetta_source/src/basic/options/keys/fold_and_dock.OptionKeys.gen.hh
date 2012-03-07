// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/fold_and_dock.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_fold_and_dock_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_fold_and_dock_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace fold_and_dock { extern BooleanOptionKey const fold_and_dock; }
namespace fold_and_dock { extern BooleanOptionKey const move_anchor_points; }
namespace fold_and_dock { extern BooleanOptionKey const set_anchor_at_closest_point; }
namespace fold_and_dock { extern BooleanOptionKey const rotate_anchor_to_x; }
namespace fold_and_dock { extern RealOptionKey const trans_mag_smooth; }
namespace fold_and_dock { extern RealOptionKey const rot_mag_smooth; }
namespace fold_and_dock { extern RealOptionKey const rb_rot_magnitude; }
namespace fold_and_dock { extern RealOptionKey const rb_trans_magnitude; }
namespace fold_and_dock { extern IntegerOptionKey const rigid_body_cycles; }
namespace fold_and_dock { extern RealOptionKey const rigid_body_frequency; }
namespace fold_and_dock { extern BooleanOptionKey const rigid_body_disable_mc; }
namespace fold_and_dock { extern RealOptionKey const slide_contact_frequency; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
