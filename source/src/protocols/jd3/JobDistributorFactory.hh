// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/JobDistributorFactory.hh
/// @brief  Definition of the JobDistributorFactory class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_JobDistributorFactory_hh
#define INCLUDED_protocols_jd3_JobDistributorFactory_hh

// Unit headers
#include <protocols/jd3/JobDistributorFactory.fwd.hh>

// Package headers
#include <protocols/jd3/JobDistributor.fwd.hh>

namespace protocols {
namespace jd3 {

class JobDistributorFactory {
public:

	static
	JobDistributorOP
	create_job_distributor();

};

} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_JobDistributorFactory_HH
