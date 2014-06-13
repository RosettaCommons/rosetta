// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/make_rot_lib/MakeRotLibMover.hh
/// @brief Forward header file for RotData
/// @author P. Douglas Renfrew ( renfrew@nyu.edu )

#ifndef INCLUDED_protocols_make_rot_lib_rotdata_fwd_hh
#define INCLUDED_protocols_make_rot_lib_rotdata_fwd_hh

#include <utility/vector1.hh>

namespace protocols {
namespace make_rot_lib {

class RotData;
typedef utility::vector1<RotData> RotDataVec;
typedef RotDataVec RotVec;

} // namespace make_rot_lib
} // namespace protocols

#endif // INCLUDED_protocols_makerotlib_rotdata_fwd_hh
