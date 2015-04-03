// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/ScreeningJobInputter.hh
/// @brief  The ScreeningJobInputter forward headers
/// @author Sam DeLuca <samuel.l.deluca@vanderbilt.edu>


#ifndef INCLUDED_protocols_jd2_ScreeningJobInputter_hh
#define INCLUDED_protocols_jd2_ScreeningJobInputter_hh

#include <protocols/jd2/ScreeningJobInputter.fwd.hh>
#include <protocols/jd2/JobInputter.hh>
#include <protocols/jd2/Job.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>


namespace protocols {
namespace jd2 {

class ScreeningJobInputter : public protocols::jd2::JobInputter
{
public:
	ScreeningJobInputter();

	virtual ~ScreeningJobInputter();

	/// @brief Fill the pose reference with the pose indicated by the job
	virtual void pose_from_job(core::pose::Pose & pose, JobOP job);

	/// @brief fill the jobs based on the specified json file
	virtual void fill_jobs(Jobs & jobs);

	/// @brief return the input source
	virtual JobInputterInputSource::Enum input_source() const;
};

}
}


#endif /* SCREENINGJOBINPUTTER_FWD_HH_ */
