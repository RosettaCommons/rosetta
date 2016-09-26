// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/JobDistributor.cc
/// @brief  JobDistributor class declaration
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_JobDistributor_FWD_HH
#define INCLUDED_protocols_jd3_JobDistributor_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace jd3 {

class JobDistributor;

typedef utility::pointer::shared_ptr< JobDistributor > JobDistributorOP;
typedef utility::pointer::shared_ptr< JobDistributor const > JobDistributorCOP;

}//jd2
}//protocols

#endif //INCLUDED_protocols_jd3_JobDistributor_HH
