// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/LoopMoverFactoryFactory.cc
/// @brief  Factory for creating LoopMovers objects
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/loops/LoopMoverFactory.hh>
#include <protocols/loops/LoopMover.hh>
#include <protocols/moves/MoverFactory.hh>

// Package Headers
#include <basic/Tracer.hh>

// Project Headers
#include <protocols/loops/Loops.hh>
#include <utility/vector0.hh>
#include <utility/exit.hh>


// C++ Headers
#include <sstream>

//Auto Headers
#include <utility/vector1.hh>

namespace protocols {
namespace loops {

using std::endl;
using std::string;
using std::pair;
using std::stringstream;
using core::pose::Pose;

static basic::Tracer tr("protocols.loops.LoopMoverFactory");

LoopMoverFactory * LoopMoverFactory::instance_( 0 );

/// @details Private constructor insures correctness of singleton.
LoopMoverFactory::LoopMoverFactory() {}

LoopMoverFactory::LoopMoverFactory(
	const LoopMoverFactory &
) {}

LoopMoverFactory::~LoopMoverFactory() {}


LoopMoverFactory *
LoopMoverFactory::get_instance()
{
	if ( instance_ == 0 ) {
		instance_ = new LoopMoverFactory;
	}
	return instance_;
}


LoopMoverOP
LoopMoverFactory::create_loop_mover(
	std::string const & type_name_in,
	Loops const & loops
) {

	std::string type_name;
	// deprecated names
	if(type_name_in == "quick_ccd"){
		type_name = "LoopMover_Perturb_QuickCCD";
	} else if(type_name_in == "sdwindow"){
		type_name = "LoopMover_SlidingWindow";
	} else if(type_name_in == "quick_ccd_moves"){
		type_name = "LoopMover_Perturb_QuickCCD_Moves";
	} else if(type_name_in == "perturb_ccd"){
		type_name = "LoopMover_Perturb_CCD";
	} else if(type_name_in == "perturb_kic"){
		type_name = "LoopMover_Perturb_KIC";
	} else {
		type_name = type_name_in;
	}

	tr.Trace << "generate LoopMover of type " << type_name << std::endl;
	LoopMoverOP loop_mover( dynamic_cast<LoopMover *>((moves::MoverFactory::get_instance()->newMover(type_name)).get()));
	if(!loop_mover){
		stringstream error_msg;
		error_msg
			<< "Attempting to create Mover "
			<< "'" << type_name << "' that is not a LoopMover." << endl
			<< "check spelling or "
			<< "register a new LoopMover in the MoverFactory" << endl;
		utility_exit_with_message(error_msg.str());
	}

	loop_mover->non_OP_loops(loops);

	return loop_mover;
}


} // namespace
} // namespace
