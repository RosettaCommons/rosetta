// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file PoseEvaluator
/// @brief PoseEvaluator
/// @detailed
///
///
/// @author Oliver Lange



// Unit Headers
#include <protocols/evaluation/PoseMetricEvaluator.hh>

// Package Headers

// Project Headers
#include <core/io/silent/SilentStruct.hh>
#include <core/pose/Pose.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <core/scoring/rms_util.hh>
// C++ headers


namespace protocols {
namespace evaluation {

using namespace core;
ConstraintEvaluator::ConstraintEvaluator( std::string tag, ConstraintSet const& )
	: SingleValuePoseEvaluator< Real >( tag )
{
}

ConstraintEvaluator::ConstraintEvaluator( std::string tag, ConstraintCOPs const& )
	: SingleValuePoseEvaluator< Real >( tag )
{
}



}
}
