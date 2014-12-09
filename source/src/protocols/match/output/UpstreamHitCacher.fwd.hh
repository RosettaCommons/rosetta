// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/output/UpstreamHitCacher.fwd.hh
/// @brief  Forward declaration for class to cache and recall the conformations of upstream hits.
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_output_UpstreamHitCacher_fwd_hh
#define INCLUDED_protocols_match_output_UpstreamHitCacher_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace match {
namespace output {

class UpstreamHitCacher;

typedef utility::pointer::shared_ptr< UpstreamHitCacher > UpstreamHitCacherOP;
typedef utility::pointer::shared_ptr< UpstreamHitCacher const > UpstreamHitCacherCOP;


}
}
}

#endif
