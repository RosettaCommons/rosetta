// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/PoseOutputter.hh
/// @brief  Definition of the %PoseOutputter class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_PoseOutputter_HH
#define INCLUDED_protocols_jd3_PoseOutputter_HH

//unit headers
#include <protocols/jd3/PoseOutputter.fwd.hh>

//package headers
#include <protocols/jd3/LarvalJob.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace jd3 {

/// @brief The %PoseOutputter
class PoseOutputter : utility::pointer::ReferenceCount
{
public:

	PoseOutputter();
	virtual ~PoseOutputter();

	virtual
	bool job_has_already_completed( LarvalJob const & job ) const = 0;

	virtual
	void mark_job_as_having_started( LarvalJob const & job ) const = 0;

	virtual
	void write_output_pose( LarvalJob const & job, core::pose::Pose const & pose ) = 0;

}; // PoseInputter

} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_PoseOutputter_HH
