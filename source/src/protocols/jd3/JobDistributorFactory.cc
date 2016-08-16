// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/JobDistributorFactory.cc
/// @brief  JobDistributorFactory class; reads the command line to create a JobDistributor
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/jd3/JobDistributorFactory.hh>

// Package headers
#include <protocols/jd3/JobDistributor.hh>

//#include <protocols/jd3/job_distributors/FileSystemJobDistributor.hh>

namespace protocols {
namespace jd3 {

JobDistributorOP
JobDistributorFactory::create_job_distributor()
{
	// TEMP!
	return JobDistributorOP( new JobDistributor );
}


} // namespace jd3
} // namespace protocols

