// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/CompletedJobOutput.hh
/// @brief  typedef for the return type of Job::run
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_CompletedJobOutput_HH
#define INCLUDED_protocols_jd3_CompletedJobOutput_HH

// Unit headers
#include <protocols/jd3/CompletedJobOutput.fwd.hh>

// Package headers
#include <protocols/jd3/Job.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace jd3 {

/// @brief Each Job will return a JobStatus and a list of JobSummary/JobResult pairs.
/// A pair's index in th job_results vector is used as its identifier for when
/// the JobResult is used as input for another job, or when it is to be output
/// by the JobQueen.
struct CompletedJobOutput
{
public:
	CompletedJobOutput() : status( jd3_job_status_failed_do_not_retry ) {}

	JobStatus status;
	utility::vector1< SingleJobOutputPair > job_results;
};

}
}

#endif
