// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/Job.hh
/// @brief  The declaration for class protocols::jd3::Job and the JobStatus enum
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_Job_FWD_HH
#define INCLUDED_protocols_jd3_Job_FWD_HH

//utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace jd3 {

class Job;

typedef utility::pointer::shared_ptr< Job > JobOP;
typedef utility::pointer::shared_ptr< Job const > JobCOP;

enum JobStatus {
	jd3_job_status_success,
	jd3_job_status_inputs_were_bad,
	jd3_job_status_failed_w_exception,
	jd3_job_status_failed_retry,
	jd3_job_status_failed_max_retries,
	jd3_job_status_failed_do_not_retry
};

} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_Job_FWD_HH
