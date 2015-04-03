// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @details
///
///
/// @author Zhe Zhang


// Unit Headers
#include <devel/replica_docking/WrapFilterAsEvaluator.hh>

// Package Headers

#include <devel/replica_docking/InteractionScoreFilter.hh>
#include <devel/replica_docking/InteractionScoreFilter.fwd.hh>
// Project Headers
#include <core/io/silent/SilentStruct.hh>
#include <core/pose/Pose.hh>

// ObjexxFCL Headers
#include <utility/tag/Tag.hh> // for test of tag
//#include <devel/replica_docking/AddEncounterConstraintMover.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/MemTracer.hh>


#include <utility/vector0.hh>
#include <utility/vector1.hh>

// C++ headers

static thread_local basic::Tracer tr( "devel.replica_docking.WrapFilterAsEvaluator" );

namespace devel {
namespace replica_docking {

using namespace core;
WrapFilterAsEvaluator::WrapFilterAsEvaluator( std::string tag ) :
	protocols::evaluation::SingleValuePoseEvaluator<core::Real>( tag )
	//	iscfilter_( scorefxn_name )
{
	irmsfilter_ = new devel::replica_docking::InteractionScoreFilter();
	tr << " Interaction Score/I_sc Evaluator " << std::endl;

}

Real WrapFilterAsEvaluator::apply( pose::Pose& pose ) const {
	return irmsfilter_->report_sm( pose );
}

bool WrapFilterAsEvaluator::applicable( pose::Pose const& pose ) const {
  return pose.is_fullatom();
}

}
}
