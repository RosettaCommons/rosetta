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
/// @details
///
///
/// @author Oliver Lange


// Unit Headers
#include <protocols/constraints_additional/ConstraintEvaluator.hh>
#include <protocols/constraints_additional/CombinedConstraintEvaluator.hh>

// Package Headers
#include <protocols/constraints_additional/MaxSeqSepConstraintSet.hh>

// Project Headers
#include <core/io/silent/SilentStruct.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <core/scoring/ScoreType.hh>


// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


// C++ headers

static THREAD_LOCAL basic::Tracer tr( "protocols.constraints_additional.ConstraintEvaluator" );

namespace protocols {
namespace constraints_additional {

using namespace core;
using namespace scoring;
using namespace constraints;


CombinedConstraintEvaluator::CombinedConstraintEvaluator( std::string tag, std::string file_name, Size combine_ratio, Size repeat )
: name_( tag )
{
	for ( Size i = 1; i <= repeat; ++i ) {
		cst_lib_.push_back( ConstraintEvaluator( tag, file_name ) );
		cst_lib_.back().set_combine_ratio( combine_ratio );
	}
}

Real CombinedConstraintEvaluator::apply( core::pose::Pose& pose_in ) const {
	Real sum( 0.0 );
	Size N( cst_lib_.size() );
	for ( Size i = 1; i <= N; ++i ) {
		sum += cst_lib_[ i ].apply(pose_in);
	}
	return sum / N;
}

std::string CombinedConstraintEvaluator::name( core::Size i ) const {
	if ( i == 1 ) { return name_; }
	if ( i == 2 ) { return name_ + "_std"; }
	runtime_assert( i <= 2 && i > 0 );
	return ""; //make compiler happy
}

void CombinedConstraintEvaluator::apply( core::pose::Pose& pose_in, std::string, core::io::silent::SilentStruct &pss ) const {
	Real sum( 0.0 );
	Real std( 0.0 );
	Size N( cst_lib_.size() );
	for ( Size i = 1; i <= N; ++i ) {
		Real val= cst_lib_[ i ].apply(pose_in);
		sum+=val;
		std+=val*val;
	}
	Real mean( sum / N );
	pss.add_energy( name( 1 ), mean );
	if ( size() == 2 ) pss.add_energy( name( 2 ), std/N-mean*mean);
}

}
}
