// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/pose_sewing/strong_types.hh
/// @author Jared Adolf-Bryfogle, jadolfbr@gmail.com

#ifndef INCLUDED_protocols_pose_sewing_StrongTypes_HH
#define INCLUDED_protocols_pose_sewing_StrongTypes_HH

//utility headers
#include <utility/strong_aliasing.hh>

namespace protocols {
namespace pose_sewing {

//Block number
using Block = utility::StrongSize< struct Block_ >;


#ifdef SERIALIZATION
SERIALIZE_STRONG_SIZE_HH( Block );
#endif

} // namespace pose_sewing
} // namespace protocols

#endif //INCLUDED_protocols_pose_sewign_StrongTypes_HH
