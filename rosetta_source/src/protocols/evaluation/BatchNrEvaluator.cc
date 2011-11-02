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
#include <protocols/evaluation/BatchNrEvaluator.hh>

// Package Headers
#include <protocols/jd2/JobDistributor.hh>

// Project Headers
#include <core/io/silent/SilentStruct.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>

// ObjexxFCL Headers
// AUTO-REMOVED #include <ObjexxFCL/string.functions.hh>

// Utility headers
// AUTO-REMOVED #include <basic/Tracer.hh>
// AUTO-REMOVED #include <core/scoring/rms_util.hh>

#include <utility/vector1.hh>

// C++ headers


namespace protocols {
namespace evaluation {

using namespace core;

core::Size
BatchNrEvaluator::apply(
 pose::Pose&
) const {
	using namespace protocols::jd2;
	return JobDistributor::get_instance()->current_batch_id();
}

void BatchEvaluator::apply( core::pose::Pose&, std::string, core::io::silent::SilentStruct &pss ) const {
	using namespace protocols::jd2;
	pss.add_string_value( name( 1 ), JobDistributor::get_instance()->get_current_batch() );
}


}
}
