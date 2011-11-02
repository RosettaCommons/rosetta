// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief demo program for implementing loop relax + FA relax
/// @author Srivatsan Raman
/// @author James Thompson
/// @author Mike Tyka
/// @author Daniel J. Mandell

// include these first for building on Visual Studio
#include <protocols/loops/Loops.hh>
#include <protocols/loops/LoopMover.fwd.hh>
#include <protocols/loops/LoopMover.hh>
#include <protocols/loops/LoopMover_SlidingWindow.hh>
#include <protocols/loops/LoopMover_QuickCCD.hh>
#include <protocols/loops/LoopMover_QuickCCD_Moves.hh>
#include <protocols/loops/LoopMover_CCD.hh>
#include <protocols/loops/LoopMover_KIC.hh>
// AUTO-REMOVED #include <protocols/loops/looprelax_protocols.hh>
// AUTO-REMOVED #include <protocols/loops/LoopBuild.hh>

#include <utility/exit.hh>

#include <string>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>



namespace protocols {
namespace loops {

static basic::Tracer TR("protocols.loops");

loops::LoopMoverOP get_loop_mover(
	std::string const & name, loops::Loops const & loops
) {
	loops::LoopMoverOP remodel_mover;
	if ( name == "quick_ccd" ) {
		remodel_mover = new loops::LoopMover_Perturb_QuickCCD( loops );
	} else if ( name == "sdwindow" ) {
		remodel_mover = new loops::LoopMover_SlidingWindow( loops );
	} else if ( name == "quick_ccd_moves" ) {
		remodel_mover = new loops::LoopMover_Perturb_QuickCCD_Moves( loops );
	} else if ( name == "perturb_ccd" ) {
		remodel_mover = new loops::LoopMover_Perturb_CCD( loops );
	} else if ( name == "perturb_kic" ) {
		remodel_mover = new loops::LoopMover_Perturb_KIC( loops );
	} else {
		std::string msg( "No mover corresponding to name " + name );
		utility_exit_with_message( msg );
	}

	return remodel_mover;
}


} // namespace loops
} // namespace protocols
