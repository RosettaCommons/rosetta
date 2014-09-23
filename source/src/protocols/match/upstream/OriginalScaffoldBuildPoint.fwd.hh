// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/OriginalScaffoldBuildPoint.fwd.hh
/// @brief  Class forward declarations for the launch point geometry on the Scaffold.
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_upstream_OriginalScaffoldBuildPoint_fwd_hh
#define INCLUDED_protocols_match_upstream_OriginalScaffoldBuildPoint_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace match {
namespace upstream {

class ProteinBackboneBuildPoint;
typedef utility::pointer::shared_ptr< ProteinBackboneBuildPoint > ProteinBackboneBuildPointOP;
typedef utility::pointer::shared_ptr< ProteinBackboneBuildPoint const > ProteinBackboneBuildPointCOP;

class OriginalBackboneBuildPoint;
typedef utility::pointer::shared_ptr< OriginalBackboneBuildPoint > OriginalBackboneBuildPointOP;
typedef utility::pointer::shared_ptr< OriginalBackboneBuildPoint const > OriginalBackboneBuildPointCOP;

}
}
}

#endif
