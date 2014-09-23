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

// Package headers
#include <protocols/loops/LoopsFileIO.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>
#include <protocols/moves/MoverFactory.hh>

// Project Headers
#include <protocols/loops/Loops.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/vector0.hh>
#include <utility/exit.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

// C++ Headers
#include <sstream>


namespace protocols {
namespace loops {

using std::endl;
using std::string;
using std::pair;
using std::stringstream;
using core::pose::Pose;

static thread_local basic::Tracer tr( "protocols.loops.LoopMoverFactory" );

#if defined MULTI_THREADED && defined CXX11
std::atomic< LoopMoverFactory * > LoopMoverFactory::instance_( 0 );
#else
LoopMoverFactory * LoopMoverFactory::instance_( 0 );
#endif

#ifdef MULTI_THREADED
#ifdef CXX11

std::mutex LoopMoverFactory::singleton_mutex_;

std::mutex & LoopMoverFactory::singleton_mutex() { return singleton_mutex_; }

#endif
#endif

/// @brief static function to get the instance of ( pointer to) this singleton class
LoopMoverFactory * LoopMoverFactory::get_instance()
{
	boost::function< LoopMoverFactory * () > creator = boost::bind( &LoopMoverFactory::create_singleton_instance );
	utility::thread::safely_create_singleton( creator, instance_ );
	return instance_;
}

LoopMoverFactory *
LoopMoverFactory::create_singleton_instance()
{
	return new LoopMoverFactory;
}

/// @details Private constructor insures correctness of singleton.
LoopMoverFactory::LoopMoverFactory() {}

LoopMoverFactory::LoopMoverFactory(
	const LoopMoverFactory &
) {}

LoopMoverFactory::~LoopMoverFactory() {}

loop_mover::LoopMoverOP
LoopMoverFactory::create_loop_mover(
	std::string const & type_name_in,
	LoopsOP const loops
)
{
	loop_mover::LoopMoverOP loop_mover = create_loop_mover( type_name_in );

	loop_mover->set_guarded_loops_not_in_charge(); // <-- tell the loop mover someone else has already resolved the loop indices.
	loop_mover->loops(loops);
	return loop_mover;
}

/// @details Set the LoopsFileData for the LoopMover leaving it in an "in
/// charge" state.  It will, upon a future call to apply, resolve the loop
/// indices, which may have been provided as PDB indices, into Pose indices.
loop_mover::LoopMoverOP
LoopMoverFactory::create_loop_mover(
	std::string const & type_name_in,
	LoopsFileData const & loops
) {
	loop_mover::LoopMoverOP loop_mover = create_loop_mover( type_name_in );
	loop_mover->loops(loops); // <-- this call leaves the loops_mover "in charge" of index resolution.
	return loop_mover;
}

loop_mover::LoopMoverOP
LoopMoverFactory::create_loop_mover(
	std::string const & type_name_in,
	GuardedLoopsFromFileOP guarded_loops
) {
	loop_mover::LoopMoverOP loop_mover = create_loop_mover( type_name_in );
	loop_mover->loops(guarded_loops); // <-- this call makes a pointer assignment.
	return loop_mover;
}


loop_mover::LoopMoverOP
LoopMoverFactory::create_loop_mover(
	std::string const & type_name_in
)
{
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
	loop_mover::LoopMoverOP loop_mover( utility::pointer::dynamic_pointer_cast< loop_mover::LoopMover > ( (moves::MoverFactory::get_instance()->newMover(type_name)) ));
	if(!loop_mover){
		stringstream error_msg;
		error_msg
			<< "Attempting to create Mover "
			<< "'" << type_name << "' that is not a LoopMover." << endl
			<< "check spelling or "
			<< "register a new LoopMover in the MoverFactory" << endl;
		utility_exit_with_message(error_msg.str());
	}
	return loop_mover;
}

} // namespace
} // namespace
