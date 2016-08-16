// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/output/MatchGrouper.fwd.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_output_MatchGrouper_fwd_hh
#define INCLUDED_protocols_match_output_MatchGrouper_fwd_hh

// Unit headers
// you cannot #include yourself #include <protocols/match/output/MatchGrouper.fwd.hh>

// Utility headers

// C++ headers

#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace match {
namespace output {

class MatchGrouper;

typedef utility::pointer::shared_ptr< MatchGrouper > MatchGrouperOP;
typedef utility::pointer::shared_ptr< MatchGrouper const > MatchGrouperCOP;


}
}
}

#endif
