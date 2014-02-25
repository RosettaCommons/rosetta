// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_ResidueSampler.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das
/// @author Parin Sripakdeevong


#ifndef INCLUDED_protocols_stepwise_PoseFilter_HH
#define INCLUDED_protocols_stepwise_PoseFilter_HH

#include <core/pose/Pose.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <protocols/stepwise/sampling/protein/PoseFilter.fwd.hh>

//Auto Headers
#include <utility/vector1.hh>


namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {

	class PoseFilter: public utility::pointer::ReferenceCount {
	public:

    PoseFilter() {}

    ~PoseFilter() {}

		virtual
		bool passes_filter( core::pose::Pose & ) = 0;

  };

} //protein
} //sampling
} //stepwise
} //protocols

#endif

