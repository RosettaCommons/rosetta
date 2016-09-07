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
#include <protocols/evaluation/PoseEvaluator.hh>

// Package Headers

// Project Headers
#include <core/io/silent/SilentStruct.hh>
#include <core/pose/Pose.hh>

// ObjexxFCL Headers

// Utility headers

#include <utility/exit.hh>
#include <utility/vector1.hh>


// C++ headers


namespace protocols {
namespace evaluation {
using namespace core;

void PoseEvaluator::apply(  io::silent::SilentStruct &pss) const {
	core::pose::Pose pose;
	pss.fill_pose( pose );
	apply( pose, pss.decoy_tag(), pss );
}

void MetaPoseEvaluator::apply( pose::Pose& pose, std::string tag, io::silent::SilentStruct &pss) const {
	for (const auto & evaluator : evaluators_) {
		evaluator->apply( pose, tag, pss );
	}
}

std::string MetaPoseEvaluator::name( core::Size ind ) const {
	runtime_assert( ind <= size() );
	Size s( 0 );
	Size last_s( 0 );
	auto it = evaluators_.begin();
	for ( ; it != evaluators_.end() && s < ind; ++it ) {
		last_s = s;
		s += (*it)->size();
	}
	runtime_assert( it != evaluators_.end() );
	return (*it)->name( ind - last_s );
}


} //evaluation
} //protocol

