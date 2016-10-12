// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_trajectory_DbTrajectoryReader_hh
#define INCLUDED_protocols_trajectory_DbTrajectoryReader_hh

// Unit Headers
#include <protocols/trajectory/DbTrajectoryReader.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>
#include <boost/noncopyable.hpp>

// C++ headers
#include <string>

namespace protocols {
namespace trajectory {

class DbTrajectoryReader
	: public utility::pointer::ReferenceCount, private boost::noncopyable {

public:

	/// @brief Default constructor.
	DbTrajectoryReader();

	/// @brief Constructor which accepts a job id.
	DbTrajectoryReader(core::Size job_id);

	/// @brief Set the job id.  This controls which records are selected.
	void set_job_id(core::Size job_id);

	/// @brief Return the number of iterations that were recorded for this job.
	core::Size get_num_iterations() const;

	/// @brief Return a vector listing every iteration recorded for this job.
	utility::vector1<core::Size> get_iterations() const;

	/// @brief Return the pose recorded for the given iteration.
	core::pose::Pose get_pose(core::Size iteration) const;

	/// @brief Return all the poses contained in this trajectory.
	utility::vector1<core::pose::Pose> get_poses() const;

private:
	utility::sql_database::sessionOP db_session_;
	core::Size job_id_;

};

} // trajectory namespace
} // protocols namespace

#endif
