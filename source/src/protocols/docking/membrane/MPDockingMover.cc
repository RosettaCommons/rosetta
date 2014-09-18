// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/MPDockingMover.fwd.hh
/// @brief      Dock two membrane proteins
/// @author     JKLeman (julia.koehler1982@gmail.com)
/// @note       Last Modified (6/24/14)

#ifndef INCLUDED_protocols_docking_membrane_MPDockingMover_cc
#define INCLUDED_protocols_docking_membrane_MPDockingMover_cc

// Unit Headers
#include <protocols/docking/membrane/MPDockingMover.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/docking/DockMCMProtocol.hh>
#include <protocols/moves/MoverContainer.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMolMover.hh>

// Package Headers
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh> 

// Utility Headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/membrane_new.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>

static thread_local basic::Tracer TR( "protocols.docking.membrane.MPDockingMover" );

namespace protocols {
namespace docking {
namespace membrane {

using namespace core;
using namespace core::pose;
using namespace protocols::moves;
		
/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default Constructor
/// @details Create a membrane pose setting the membrane center
/// at center=(0, 0, 0), normal=(0, 0, 1) and loads in spans
/// and lips from the command line interface.
MPDockingMover::MPDockingMover() :
	protocols::moves::Mover(),
	center_(0, 0, 0),
	normal_(0, 0, 1)
{}

/// @brief Copy Constructor
/// @details Create a deep copy of this mover
MPDockingMover::MPDockingMover( MPDockingMover const & src ) :
	protocols::moves::Mover( src ), 
	center_( src.center_ ),
	normal_( src.normal_ )
{}

/// @brief Destructor
MPDockingMover::~MPDockingMover() {}

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
MPDockingMover::clone() const {
	return ( new MPDockingMover( *this ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
MPDockingMover::fresh_instance() const {
	return new MPDockingMover();
}

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this Mover (MPDockingMover)
std::string
MPDockingMover::get_name() const {
	return "MPDockingMover";
}


/// @brief Add Membrane Components to Pose
/// @details Add membrane components to pose which includes
///	spanning topology, lips info, embeddings, and a membrane
/// virtual residue describing the membrane position
void
MPDockingMover::apply( Pose & pose ) {
	
	using namespace core::conformation::membrane;

	TR << "calling setup" << std::endl;
	setup();

	// assuming that protein 1 is fixed in the membrane!!!
	// add membrane VRT, call AddMembraneMover
	TR << "adding MEM" << std::endl;
	add_membrane_mover_->apply( pose );

	// foldtree
	TR << "creating foldtree" << std::endl;
	core::kinematics::FoldTree foldtree = pose.fold_tree();
	foldtree.reorder( 81 );
	pose.fold_tree( foldtree );

	// add a jump from protein 1 (fixed) to protein 2 (flexible)
//	foldtree.add_edge( protein2_start, protein2_end, number );
	
	// check foldtree
	TR << foldtree << std::endl;;
	
	// create MC-object
	TR << "create MC object" << std::endl;
	protocols::moves::MonteCarloOP montecarlo = new protocols::moves::MonteCarlo( pose, *scorefunction_, kT_);
	
	// attach Pymol observer
	TR << "attach Pymol observer" << std::endl;
	protocols::moves::AddPyMolObserver( pose );

	// make a move, rigidbody mover
	TR << "calling dock_MCM_protocol" << std::endl;
	dock_mcm_protocol_->apply( pose );

	// score
	TR << "calling Metropolis" << std::endl;
	montecarlo->boltzmann( pose );
		
}

////////////////////////////////////////////////////////////////////////////////
void MPDockingMover::setup(){
	add_membrane_mover_ = new protocols::membrane::AddMembraneMover();
	scorefunction_ = core::scoring::getScoreFunction();
	dock_mcm_protocol_ = new docking::DockMCMProtocol( 1, scorefunction_, scorefunction_ );
	random_mover_ = new RandomMover();
//	scorefunction_ = core::scoring::createScoreFunction( "cen_membrane_2014.wts" );
	kT_ = 1;

}

} // membrane
} // docking
} // protocols

#endif // INCLUDED_protocols_docking_membrane_MPDockingMover_cc
