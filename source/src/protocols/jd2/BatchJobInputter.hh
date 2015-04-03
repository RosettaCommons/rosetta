// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/BatchJobInputter.hh
/// @brief  header file for BatchJobInputter class
/// @author Oliver Lange


#ifndef INCLUDED_protocols_jd2_BatchJobInputter_hh
#define INCLUDED_protocols_jd2_BatchJobInputter_hh

//unit headers
#include <protocols/jd2/JobInputter.hh>
#include <protocols/jd2/BatchJobInputter.fwd.hh>
#include <protocols/jd2/Job.hh>

//project headers
#include <core/pose/Pose.fwd.hh>
#include <utility/options/OptionCollection.hh>

#include <utility/vector1.hh>


//utility headers

namespace protocols {
namespace jd2 {

/// @details This is the simplest implementation of JobInputter, which reads from -s/-l and Batch files.
class BatchJobInputter : public protocols::jd2::JobInputter
{
public:
	static std::string const BOGUS_BATCH_ID;
	BatchJobInputter( std::string batch );

	virtual ~BatchJobInputter();

	/// @brief this function is responsible for filling the pose reference with
	/// the pose indicated by the job.  The Job object (within its InnerJob)
	/// contains a PoseCOP.  This function needs to either fill the pose
	/// reference from the InnerJob or, on first demand of a pose from that
	/// InnerJob, instantiate the pose, hand off a COP to the InnerJob, and fill
	/// the reference.  This implementation uses pose_from_pdb
	virtual void pose_from_job( core::pose::Pose & pose, JobOP job ) {
		check_batch();
		this_batch_job_inputter_->pose_from_job( pose, job);
	}

	/// @brief this function determines what jobs exist from -in::file::silent and
	/// -in::file::tags.
	virtual void fill_jobs( Jobs & jobs ) {
		check_batch();
		this_batch_job_inputter_->fill_jobs( jobs );
	}

	/// @brief Return the type of input source that the BatchJobInputter is currently
	///  using.
	/// @return The input source for the current batch.
	virtual JobInputterInputSource::Enum input_source() const;

private:
	void read_batch();
	void check_batch();
	std::string current_batch_;
	JobInputterOP this_batch_job_inputter_;
	utility::options::OptionCollection const vanilla_options_; //options before batch-options were added
}; // BatchJobInputter

} // namespace jd2
} // namespace protocols

#endif //INCLUDED_protocols_jd2_BatchJobInputter_HH
