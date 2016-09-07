// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @details
///   Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange


// Unit Headers
#include <protocols/evaluation/TimeEvaluator.hh>

// Package Headers
#include <protocols/evaluation/PoseEvaluator.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


// ObjexxFCL Headers

// Utility headers

//// C++ headers

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif

namespace protocols {
namespace evaluation {


TimeEvaluator::TimeEvaluator( std::string const& tag )
: evaluation::SingleValuePoseEvaluator< core::Real > ("time"+tag)
{
	reset();
}

void
TimeEvaluator::reset() {
	start_time_ = time(nullptr);
}


core::Real
TimeEvaluator::apply( core::pose::Pose& ) const
{
	return time(nullptr) - start_time_;
}

}
}
