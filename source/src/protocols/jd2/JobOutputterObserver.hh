// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/JobOutputterObserver.hh
/// @brief  classes that want to add score_value pairs to poses before they are output should derive from this Observer class
/// @author Oliver Lange


#ifndef INCLUDED_protocols_jd2_JobOutputterObserver_hh
#define INCLUDED_protocols_jd2_JobOutputterObserver_hh
//project headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>


namespace protocols {
namespace jd2 {

class JobOutputterObserver {
public:
	virtual
	void add_values_to_job( core::pose::Pose const& pose, protocols::jd2::JobOP ) const = 0;
};

#include <utility/pointer/access_ptr.hh>
typedef utility::pointer::access_ptr< JobOutputterObserver const > JobOutputterObserverAP;

} // namespace jd2
} // namespace protocols


#endif //INCLUDED_protocols_jd2_JobOutputterObserver_HH

