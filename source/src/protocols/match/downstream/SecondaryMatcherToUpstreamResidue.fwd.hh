// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/downstream/DownstreamAlgorithm.fwd.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_downstream_SecondaryMatcherToUpstreamResidue_fwd_hh
#define INCLUDED_protocols_match_downstream_SecondaryMatcherToUpstreamResidue_fwd_hh

// Unit headers
// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace match {
namespace downstream {

class SecondaryMatcherToUpstreamResidue;
typedef utility::pointer::shared_ptr< SecondaryMatcherToUpstreamResidue > SecondaryMatcherToUpstreamResidueOP;
typedef utility::pointer::shared_ptr< SecondaryMatcherToUpstreamResidue const > SecondaryMatcherToUpstreamResidueCOP;

class TargetRotamerCoords;
typedef utility::pointer::shared_ptr< TargetRotamerCoords > TargetRotamerCoordsOP;
typedef utility::pointer::shared_ptr< TargetRotamerCoords const > TargetRotamerCoordsCOP;


class SecondaryMatchUpstreamResProcessor;
typedef utility::pointer::shared_ptr< SecondaryMatchUpstreamResProcessor > SecondaryMatchUpstreamResProcessorOP;
typedef utility::pointer::shared_ptr< SecondaryMatchUpstreamResProcessor const > SecondaryMatchUpstreamResProcessorCOP;

}
}
}

#endif
