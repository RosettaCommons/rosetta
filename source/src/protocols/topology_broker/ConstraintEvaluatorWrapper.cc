// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file TopologyBroker
/// @brief  top-class (Organizer) of the TopologyBroker mechanism
/// @details responsibilities:
/// @author Oliver Lange

// Unit Headers
#include <protocols/topology_broker/ConstraintEvaluatorWrapper.hh>

// Package Headers

// Project Headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
// ObjexxFCL Headers

// Utility headers
#include <utility>
#include <utility/excn/Exceptions.hh>

#include <basic/Tracer.hh>
//// C++ headers

#include <core/id/SequenceMapping.hh>
#include <utility/vector1.hh>


// option key includes

static basic::Tracer tr( "protocols.evaluator", basic::t_info );

namespace protocols {
namespace topology_broker {


using namespace core;
using namespace scoring::constraints;
using namespace scoring;
ConstraintEvaluatorWrapper::ConstraintEvaluatorWrapper( std::string  name, ConstraintClaimerCOP claimer ) :
	name_(std::move( name )),
	claimer_(std::move( claimer ))
{}

Real ConstraintEvaluatorWrapper::apply( core::pose::Pose& pose_in ) const {
	pose::Pose pose( pose_in );
	claimer_->add_constraints( pose );

	ScoreFunction scfxn;
	scfxn.set_weight( scoring::atom_pair_constraint, 1.0 );
	core::Real score( scfxn( pose ) );
	return score;
}

void ConstraintEvaluatorWrapper::apply( core::pose::Pose& pose_in, std::string, core::io::silent::SilentStruct &pss ) const {
	pose::Pose pose( pose_in );
	claimer_->add_constraints( pose );

	ScoreFunction scfxn;
	scfxn.set_weight( scoring::atom_pair_constraint, 1.0 );
	core::Real score( scfxn( pose ) );

	pss.add_energy( name( 1 ), score );
}

std::string ConstraintEvaluatorWrapper::name( core::Size i ) const {
	if ( i == 1 ) { return name_; }
	runtime_assert( i <= 1 && i > 0 );
	return ""; //make compiler happy
}

} //topology_broker
} //protocols
