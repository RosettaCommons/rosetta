// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/job_distributors/MultiThreadedJobDistributor.fwd.hh
/// @brief  MultiThreadedJobDistributor class declaration
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_job_distributors_MultiThreadedJobDistributor_FWD_HH
#define INCLUDED_protocols_jd3_job_distributors_MultiThreadedJobDistributor_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace jd3 {
namespace job_distributors {

class MultiThreadedJobDistributor;
class JobRunner;

typedef utility::pointer::shared_ptr< MultiThreadedJobDistributor > MultiThreadedJobDistributorOP;
typedef utility::pointer::shared_ptr< MultiThreadedJobDistributor const > MultiThreadedJobDistributorCOP;

typedef utility::pointer::shared_ptr< JobRunner > JobRunnerOP;
typedef utility::pointer::shared_ptr< JobRunner const > JobRunnerCOP;

} // job_distributors
} // jd3
} // protocols

#endif //INCLUDED_protocols_jd3_job_distributors_MultiThreadedJobDistributor_HH
