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

// Package headers
#include <protocols/jd3/JobSummary.fwd.hh>
#include <protocols/jd3/JobResult.fwd.hh>

// C++ headers
#include <utility>

namespace protocols {
namespace jd3 {

typedef std::pair< JobSummaryOP, JobResultOP > CompletedJobOutput;

}
}

#endif
