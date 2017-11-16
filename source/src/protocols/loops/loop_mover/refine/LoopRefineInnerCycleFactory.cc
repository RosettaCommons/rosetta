// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/loop_mover/refine/LoopRefineInnerCycleFactory.cc
/// @brief  Factory for creating LoopRefineInnerCycle objects
/// @author Brian D. Weitzner ( brian.weitzner@gmail.com )

// Unit Headers
#include <protocols/loops/loop_mover/refine/LoopRefineInnerCycleFactory.hh>

// Package headers
#include <protocols/loops/loop_mover/refine/LoopRefineInnerCycle.hh>
#include <protocols/loops/loop_mover/refine/LoopRefineInnerCycleContainer.hh>
#include <protocols/moves/MoverFactory.hh>

// Project Headers
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/moves/MonteCarlo.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/vector1.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

// C++ Headers
#include <sstream>

namespace protocols {
namespace loops {
namespace loop_mover {
namespace refine {

using std::string;
using utility::tools::make_vector1;
using utility::vector1;

static basic::Tracer TR( "protocols.loops.loop_mover.refine.LoopRefineInnerCycleFactory" );

void LoopRefineInnerCycleFactory::setup_known_types()
{
	// It's ok to skip checking if the vector has been set up already because this class is a singleton.
	loop_refine_inner_cycle_name_to_string_.resize( number_of_loop_refine_inner_cycle_names );

	// Individual LoopRefineInnerCycles
	loop_refine_inner_cycle_name_to_string_[ IC_SmallMinCCDTrial ] = make_vector1< string >( "SmallMinCCDTrial" );
	loop_refine_inner_cycle_name_to_string_[ IC_ShearMinCCDTrial ] = make_vector1< string >( "ShearMinCCDTrial" );
	loop_refine_inner_cycle_name_to_string_[ IC_RepackTrial ] = make_vector1< string >( "RepackTrial" );

	// Pre-made algorithms
	loop_refine_inner_cycle_name_to_string_[ IC_RefineCCDStandard ] = make_vector1< string >(
		"SmallMinCCDTrial", "ShearMinCCDTrial", "RepackTrial" );
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// BOILER PLATE CODE //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////


/// @details Private constructor insures correctness of singleton.
LoopRefineInnerCycleFactory::LoopRefineInnerCycleFactory()
{
	setup_known_types();
}

LoopRefineInnerCycleFactory::~LoopRefineInnerCycleFactory() {}


///////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// END OF BOILER PLATE CODE //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

LoopRefineInnerCycleOP LoopRefineInnerCycleFactory::create_inner_cycle(
	LoopRefineInnerCycleName type_name,
	LoopMover_Refine_CCDAP loop_mover,
	moves::MonteCarloOP mc,
	core::scoring::ScoreFunctionOP scorefxn,
	core::pack::task::TaskFactoryOP tf
) const {

	vector1< string > movers_to_make = loop_refine_inner_cycle_name_to_string_[ type_name ];

	LoopRefineInnerCycleOP inner_cycle;
	if ( movers_to_make.size() == 1 ) {
		inner_cycle = make_inner_cycle_from_string_name( movers_to_make[ 1 ] );
	} else {
		LoopRefineInnerCycleContainerOP tmp_inner_cycle( new LoopRefineInnerCycleContainer );

		for ( vector1< string >::const_iterator it = movers_to_make.begin(); it != movers_to_make.end(); ++it ) {
			tmp_inner_cycle->add_inner_cycle_step( make_inner_cycle_from_string_name( *it ) );
		}
		inner_cycle = tmp_inner_cycle;
	}

	inner_cycle->set_loop_mover( loop_mover );
	inner_cycle->set_mc( mc );
	inner_cycle->set_scorefxn( scorefxn );
	inner_cycle->set_task_factory( tf );

	return inner_cycle;
}


LoopRefineInnerCycleOP LoopRefineInnerCycleFactory::make_inner_cycle_from_string_name( std::string const & name ) const
{
	TR.Trace << "generate LoopRefineInnerCycle of type " << name << std::endl;
	LoopRefineInnerCycleOP inner_cycle = utility::pointer::dynamic_pointer_cast< protocols::loops::loop_mover::refine::LoopRefineInnerCycle > ( ( moves::MoverFactory::get_instance()->newMover( name ) ) );

	if ( ! inner_cycle ) {
		using utility::excn::EXCN_Msg_Exception;
		throw EXCN_Msg_Exception( "Attempting to create Mover '" + name + \
			"' that cannot be casted to a LoopRefineInnerCycle.  Check your spelling and/or confirm this mover has been " +\
			"registered to the MoverFactory." );
	}
	return inner_cycle;
}

} // namespace refine
} // namespace loop_mover
} // namespace loops
} // namespace protocols
