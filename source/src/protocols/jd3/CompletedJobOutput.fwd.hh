// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/CompletedJobOutput.fwd.hh
/// @brief  Forward declaration for the stuct return type of Job::run
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_CompletedJobOutput_FWD_HH
#define INCLUDED_protocols_jd3_CompletedJobOutput_FWD_HH

// Package headers
#include <protocols/jd3/JobSummary.fwd.hh>
#include <protocols/jd3/JobResult.fwd.hh>

// Project headers
#include <core/types.hh>

// C++ headers
#include <utility>

namespace protocols {
namespace jd3 {

typedef std::pair< JobSummaryOP, JobResultOP > SingleJobOutputPair;

/// The first entry in this pair is the Job id; the second entry is
/// the index for a particular job-output-pair in the CompletedJobOutput
/// vector. E.g., the 5th JobSummary/JobResult pair for Job 132 would
/// be identified by the pair { 132, 5 }. A JobResultID is exactly
/// the same thing, but is given another name for when the JobQueen
/// is trying to indicate which JobResult of all the JobSummary/JobResult pairs
/// that she is interested in.
typedef std::pair< core::Size, core::Size > JobResultID;

struct CompletedJobOutput;

}
}

#endif
