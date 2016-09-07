// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/LargeNstructJobInputter.hh
/// @brief  Header file for a JobInputter for cases where the number of jobs is large enough to fill up memory.
/// @author Vikram K. Mulligan, Baker Laboratory (vmullig@uw.edu)


#ifndef INCLUDED_protocols_jd2_LargeNstructJobInputter_hh
#define INCLUDED_protocols_jd2_LargeNstructJobInputter_hh

//unit headers
#include <protocols/jd2/JobInputter.hh>
#include <protocols/jd2/LargeNstructJobInputter.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/jd2/JobsContainer.hh>

//project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


//utility headers

namespace protocols {
namespace jd2 {

/// @details This is an implementation of JobInputter for cases where it's not possible to list all jobs at once.
class LargeNstructJobInputter : public protocols::jd2::JobInputter
{
public:
	LargeNstructJobInputter();

	/// @brief this function is responsible for filling the pose reference with
	/// the pose indicated by the job.  The Job object (within its InnerJob)
	/// contains a PoseCOP.  This function needs to either fill the pose
	/// reference from the InnerJob or, on first demand of a pose from that
	/// InnerJob, instantiate the pose, hand off a COP to the InnerJob, and fill
	/// the reference.
	void pose_from_job( core::pose::Pose & pose, protocols::jd2::JobOP job ) override;

	/// @brief This function determines what jobs exist.
	/// @details  Unlike the GenericJobInputter, this version only ever has a subset of total jobs in
	/// memory at any given time.
	void fill_jobs( protocols::jd2::JobsContainer & jobs ) override;

	/// @brief This function is only called by certain JobInputters to update the jobs list after it has already been created.
	/// @details An example case would be the LargeNstructJobInputter, which uses this function to load additional jobs after
	/// the first N have started to come back.
	void update_jobs_list( JobsContainerOP jobs ) override;

	/// @brief Return the type of input source that the LargeNstructJobInputter is currently
	///  using.
	/// @return Always <em>POSE</em>.
	protocols::jd2::JobInputterInputSource::Enum input_source() const override;

	/// @brief Does this type of JobInputter update the jobs list?
	/// @details False by default.  Override this function in derived classes to make it true.
	/// The LargeNstructJobInputter overrides this, and returns true.
	bool updates_jobs_list() const override { return true; }


private:
	/// @brief Private function to add N jobs to the list of jobs.
	///
	void populate_next_n_jobs(
		protocols::jd2::JobsContainer & jobs,
		core::Size const first_job_index,
		core::Size const number_of_jobs_to_add,
		core::Size const total_jobs
	);
}; // LargeNstructJobInputter class

} // namespace jd2
} // namespace protocols

#endif //INCLUDED_protocols_jd2_LargeNstructJobInputter_hh
