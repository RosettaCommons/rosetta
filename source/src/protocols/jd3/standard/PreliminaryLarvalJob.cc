// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/jd3/standard/PreliminaryLarvalJob.cc
/// @brief A simple class for that stores input information for each job defintion/input nstruct as a PreliminaryLarvalJob.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#include <protocols/jd3/standard/PreliminaryLarvalJob.hh>
#include <protocols/jd3/InnerLarvalJob.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.jd3.standard.PreliminaryLarvalJob" );


namespace protocols {
namespace jd3 {
namespace standard {

PreliminaryLarvalJob::PreliminaryLarvalJob() = default;
PreliminaryLarvalJob::~PreliminaryLarvalJob() = default;
PreliminaryLarvalJob::PreliminaryLarvalJob( PreliminaryLarvalJob const & /*src*/ ) = default;

PreliminaryLarvalJob &
PreliminaryLarvalJob::operator = ( PreliminaryLarvalJob const & rhs )
{
	if ( this != &rhs ) {
		inner_job = rhs.inner_job;
		job_tag   = rhs.job_tag;
		job_options = rhs.job_options;
		pose_inputter = rhs.pose_inputter;
	}
	return *this;
}



} //protocols
} //jd3
} //standard






