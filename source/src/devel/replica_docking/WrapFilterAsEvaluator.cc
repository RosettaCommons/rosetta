// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @details
///
///
/// @author Zhe Zhang


// Unit Headers
#include <devel/replica_docking/WrapFilterAsEvaluator.hh>

// Package Headers

// Project Headers
#include <core/pose/Pose.fwd.hh>

// ObjexxFCL Headers
//#include <devel/replica_docking/AddEncounterConstraintMover.hh>

// Utility headers
#include <basic/Tracer.hh>



#include <protocols/filters/Filter.hh> // AUTO IWYU For Filter

// C++ headers

static basic::Tracer tr( "devel.replica_docking.WrapFilterAsEvaluator" );

namespace devel {
namespace replica_docking {

using namespace core;
WrapFilterAsEvaluator::WrapFilterAsEvaluator( protocols::filters::FilterOP filter, std::string tag ) :
	protocols::evaluation::SingleValuePoseEvaluator<core::Real>( tag )
	// iscfilter_( scorefxn_name )
{
	// iscfilter_ = new devel::replica_docking::InteractionScoreFilter( scorefxn_name );
	filter_ = filter;
	tr << "Evaluator " << tag << std::endl;

}

Real WrapFilterAsEvaluator::apply( pose::Pose& pose ) const {
	return filter_->report_sm( pose );
}

bool WrapFilterAsEvaluator::applicable( pose::Pose const& /*pose*/ ) const {
	return true; //pose.is_fullatom();
}

}
}
