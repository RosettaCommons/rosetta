// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/loops/LoopMoverFactory.hh
/// @brief Factory for creating LoopMover objects
/// @author Matthew O'Meara (mattjomeara@gmail.com)


#ifndef INCLUDED_protocols_loops_LoopMoverFactory_hh
#define INCLUDED_protocols_loops_LoopMoverFactory_hh

// Unit Headers
#include <protocols/loops/LoopMoverFactory.fwd.hh>

// Project Headers
#include <protocols/loops/loop_mover/LoopMover.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>

// Platform Headers
#include <core/pose/Pose.fwd.hh>
#include <utility/factory/WidgetRegistrator.hh>

// C++ Headers
#include <map>

#include <utility/vector1.hh>


namespace protocols {
namespace loops {

/// Create LoopMover Reporters
class LoopMoverFactory {

	// Private constructor to make it singleton managed
	LoopMoverFactory();
	LoopMoverFactory(const LoopMoverFactory & src); // unimplemented

	LoopMoverFactory const &
	operator=( LoopMoverFactory const & ); // unimplemented

public:

	// Warning this is not called because of the singleton pattern
	virtual ~LoopMoverFactory();

	static LoopMoverFactory * get_instance();

	loop_mover::LoopMoverOP create_loop_mover(std::string const & type_name, LoopsOP const  loops );

private:

	static LoopMoverFactory * instance_;

};

} // namespace
} // namespace

#endif
