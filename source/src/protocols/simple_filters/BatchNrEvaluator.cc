// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file PoseEvaluator
/// @brief PoseEvaluator
/// @details
///
///
/// @author Oliver Lange


// Unit Headers
#include <protocols/simple_filters/BatchNrEvaluator.hh>

// Package Headers
#include <protocols/jd2/util.hh>

// Project Headers
#include <core/io/silent/SilentStruct.hh>

// ObjexxFCL Headers

// Utility headers

#include <utility/vector1.hh>

// C++ headers


namespace protocols {
namespace simple_filters {

using namespace core;

core::Size
BatchNrEvaluator::apply(
	pose::Pose&
) const {
	return protocols::jd2::current_batch_id();
}

void BatchEvaluator::apply( core::pose::Pose&, std::string, core::io::silent::SilentStruct &pss ) const {
	pss.add_string_value( name( 1 ), protocols::jd2::current_batch() );
}


}
}
