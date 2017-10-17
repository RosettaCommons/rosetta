// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/pose_outputters/SecondaryPoseOutputter.cc
/// @brief  Definition of the %SecondaryPoseOutputter class's methods
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

//unit headers
#include <protocols/jd3/pose_outputters/SecondaryPoseOutputter.hh>

namespace protocols {
namespace jd3 {
namespace pose_outputters {

SecondaryPoseOutputter::SecondaryPoseOutputter() {}
SecondaryPoseOutputter::~SecondaryPoseOutputter() {}

void SecondaryPoseOutputter::determine_job_tag(
	utility::tag::TagCOP,
	utility::options::OptionCollection const &,
	InnerLarvalJob &
) const
{}


bool SecondaryPoseOutputter::job_has_already_completed( LarvalJob const &, utility::options::OptionCollection const & ) const
{
	return false;
}

void SecondaryPoseOutputter::mark_job_as_having_started( LarvalJob const &, utility::options::OptionCollection const & ) const
{
}


} // namespace pose_outputters
} // namespace jd3
} // namespace protocols
