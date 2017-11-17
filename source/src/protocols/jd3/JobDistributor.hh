// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/JobDistributor.cc
/// @brief  JobDistributor class definition
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_JobDistributor_HH
#define INCLUDED_protocols_jd3_JobDistributor_HH

// Unit headers
#include <protocols/jd3/JobDistributor.fwd.hh>

// Package headers
#include <protocols/jd3/JobQueen.fwd.hh>
#include <protocols/jd3/Job.fwd.hh>
#include <protocols/jd3/JobResult.fwd.hh>
#include <protocols/jd3/LarvalJob.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace jd3 {

/// @brief OK -- I don't know what the division of labor between a JobDistributor base class
/// and a JobDistributor subclass, so I'm just going to start writing a simple JobDistributor
/// and then follow it by writing an MPI job distributor and then I'll see what shakes out.
class JobDistributor : public utility::pointer::ReferenceCount {
public:

	JobDistributor();
	~JobDistributor() override;

	/// @brief The main method for executing a protocol.
	virtual
	void
	go( JobQueenOP queen ) = 0;


};

}//jd2
}//protocols

#endif //INCLUDED_protocols_jd3_JobDistributor_HH
