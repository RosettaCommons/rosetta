// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/pose_sewing/strong_types.cc
/// @author Jared Adolf-Bryfogle, jadolfbr@gmail.com

#ifdef SERIALIZATION

//utility headers
#include <protocols/pose_sewing/strong_types.hh>
#include <utility/strong_aliasing.hh>

#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#include <utility/serialization/serialization.hh>

namespace protocols {
namespace pose_sewing {

SERIALIZE_STRONG_SIZE_CC( Block );

} // namespace pose_sewing
} // namespace protocols

#endif
