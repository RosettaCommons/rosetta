// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/jd3/standard/PreliminaryLarvalJobTracker.fwd.hh
/// @brief A class that tracks PreliminaryLarvalJobs in the SJQ.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_jd3_standard_PreliminaryLarvalJobTracker_fwd_hh
#define INCLUDED_protocols_jd3_standard_PreliminaryLarvalJobTracker_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace jd3 {
namespace standard {

class PreliminaryLarvalJobTracker;

typedef utility::pointer::shared_ptr< PreliminaryLarvalJobTracker > PreliminaryLarvalJobTrackerOP;
typedef utility::pointer::shared_ptr< PreliminaryLarvalJobTracker const > PreliminaryLarvalJobTrackerCOP;

} //protocols
} //jd3
} //standard

#endif //INCLUDED_protocols_jd3_standard_PreliminaryLarvalJobTracker_fwd_hh
