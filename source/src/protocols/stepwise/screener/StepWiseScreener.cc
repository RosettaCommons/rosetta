// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/StepWiseScreener.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/StepWiseScreener.hh>
#include <protocols/moves/CompositionMover.hh>
#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "protocols.stepwise.screener.StepWiseScreener" );

namespace protocols {
namespace stepwise {
namespace screener {

	//Constructor
	StepWiseScreener::StepWiseScreener():
		utility::pointer::ReferenceCount(),
		count_( 0 ),
		ok_to_increment_( true ) // silly hack.
	{}

	//Destructor
	StepWiseScreener::~StepWiseScreener()
	{}

	void
	StepWiseScreener::add_mover( moves::CompositionMoverOP update_mover, moves::CompositionMoverOP restore_mover ){
		update_mover->add_mover( 0 );
		restore_mover->add_mover( 0 );
	}

	void
	StepWiseScreener::increment_count(){
		if ( ok_to_increment_ ) count_++;
		ok_to_increment_ = false; // needs to be manually reset by SampleAndScreen. Useful if you don't want to increment in some inner modeler loops.
	}

} //screener
} //stepwise
} //protocols
