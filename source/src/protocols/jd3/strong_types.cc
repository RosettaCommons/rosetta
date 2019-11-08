// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/strong_types.cc
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifdef SERIALIZATION

//utility headers
#include <protocols/jd3/strong_types.hh>
#include <utility/strong_aliasing.hh>

#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#include <utility/serialization/serialization.hh>

namespace protocols {
namespace jd3 {

SERIALIZE_STRONG_SIZE_CC( JobDAGNodeID );
SERIALIZE_STRONG_SIZE_CC( PrelimJobNodeID );
SERIALIZE_STRONG_SIZE_CC( GlobalJobID );
SERIALIZE_STRONG_SIZE_CC( LocalJobID );
SERIALIZE_STRONG_SIZE_CC( NStructIndex );
SERIALIZE_STRONG_SIZE_CC( ResultIndex );

} // namespace jd3
} // namespace protocols

#endif
