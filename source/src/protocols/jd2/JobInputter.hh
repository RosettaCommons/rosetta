// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/JobInputter.hh
/// @brief  header file for JobInputter class, part of August 2008 job distributor as planned at RosettaCon08
/// @author Steven Lewis smlewi@gmail.com


#ifndef INCLUDED_protocols_jd2_JobInputter_hh
#define INCLUDED_protocols_jd2_JobInputter_hh

//unit headers
#include <core/types.hh>
#include <protocols/jd2/JobInputter.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/jd2/JobsContainer.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace jd2 {


/// @brief The %JobInputter class is responsible for
/// - determining what jobs exist, and
/// - taking a job object and returning the starting pose associated with that job.
///
/// @details NOTE: your JobInputter should order Job objects in the Jobs vector to
/// have as few "transitions" between inputs as possible (group all Jobs of the
/// same input next to each other).  This improves efficiency of the "FAIL_BAD_INPUT"
/// functionality.  I said "should", not "must" on purpose.
class JobInputter : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< JobInputter >
{
public:

	~JobInputter() override;

	/// @brief Get self const owning pointers.
	///
	inline JobInputterCOP get_self_ptr() const { return shared_from_this(); }

	/// @brief Get self owning pointers.
	///
	inline JobInputterOP get_self_ptr() { return shared_from_this(); }

	/// @brief Get self const weak pointers.
	///
	inline JobInputterCAP get_self_weak_ptr() const { return JobInputterCAP( shared_from_this() ); }

	/// @brief Get self weak pointers.
	///
	inline JobInputterAP get_self_weak_ptr() { return JobInputterAP( shared_from_this() ); }


	/// @brief this function is responsible for returning the starting pose to use for
	/// the input job.  The Job object (within its InnerJob) contains a PoseCOP.
	/// This function needs to either create a new Pose from the InnerJob Pose or, on
	/// first demand of a pose from that InnerJob, instantiate the pose, hand off a COP
	/// to the InnerJob, and fill the reference.
	virtual void pose_from_job( core::pose::Pose & pose, JobOP job ) = 0;

	/// @brief this function determines what jobs exist.  This function neither knows nor
	/// cares what jobs are already complete on disk/memory - it just figures out what
	/// ones should exist given the input.  NOTE: your JobInputter should order Job objects
	/// in the Jobs vector to have as few "transitions" between inputs as possible (group
	/// all Jobs of the same input next to each other).  This improves efficiency of the
	/// "FAIL_BAD_INPUT" functionality.  Note I (Brave Sir Steven) said "should", not "must".
	virtual void fill_jobs( JobsContainer & jobs ) = 0;

	/// @brief This function is only called by certain JobInputters to update the jobs list after it has already been created.
	/// @details An example case would be the LargeNstructJobInputter, which uses this function to load additional jobs after the first N have started to come back.
	virtual void update_jobs_list( JobsContainerOP jobs );

	/// @brief return the type of input source that the JobInputter is currently
	/// using
	virtual JobInputterInputSource::Enum input_source() const = 0;

	/// @brief call this with input_source() to get the input source of a
	/// particular job inputter
	static
	std::string
	job_inputter_input_source_to_string(
		JobInputterInputSource::Enum source
	);

	/// @brief Does this type of JobInputter update the jobs list?
	/// @details False by default.  Override this function in derived classes to make it true.
	virtual bool updates_jobs_list() const { return false; }

protected:

	/// @brief this function modifies the InnerJob's pose.  Access to that pose is via friendship.
	void load_pose_into_job( core::pose::Pose const & pose, JobOP job );

	/// @brief this function modifies the InnerJob's pose.  Access to that pose is via friendship.
	void load_pose_into_job( core::pose::PoseCOP pose, JobOP job );

	virtual core::Size get_nstruct() const;

}; // JobInputter


} // namespace jd2
} // namespace protocols

#endif //INCLUDED_protocols_jd2_JobInputter_HH
