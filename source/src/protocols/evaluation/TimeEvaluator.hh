// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @details
///   Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange


#ifndef INCLUDED_protocols_evaluation_TimeEvaluator_hh
#define INCLUDED_protocols_evaluation_TimeEvaluator_hh


// Unit Headers
#include <protocols/evaluation/TimeEvaluator.fwd.hh>

// Package Headers
#include <protocols/evaluation/PoseEvaluator.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


// ObjexxFCL Headers

// Utility headers


//// C++ headers

namespace protocols {
namespace evaluation {


class TimeEvaluator : public evaluation::SingleValuePoseEvaluator< core::Real > {
public:
	TimeEvaluator( std::string const& tag = "");

	/// @brief reset start_time
	void reset();

	/// @brief evaluate pose
	virtual core::Real apply( core::pose::Pose& ) const;
private:
	time_t start_time_;
};


}
}

#endif
