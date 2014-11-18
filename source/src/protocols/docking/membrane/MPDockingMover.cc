// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/MPDockingMover.cc
/// @brief      Dock two membrane proteins
/// @author     JKLeman (julia.koehler1982@gmail.com)
/// @note       Last Modified (6/24/14)

#ifndef INCLUDED_protocols_docking_membrane_MPDockingMover_cc
#define INCLUDED_protocols_docking_membrane_MPDockingMover_cc

// Unit Headers
#include <protocols/docking/membrane/MPDockingMover.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>

#include <protocols/membrane/AddMembraneMover.hh>

#include <protocols/docking/DockingProtocol.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMolMover.hh>

// Package Headers
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
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
using namespace core::scoring;
using namespace protocols::membrane;
using namespace protocols::moves;
using namespace protocols::docking;
	
/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default Constructor
/// @details Docks two proteins with default normal=(0,0,1) and center=(0,0,0)
MPDockingMover::MPDockingMover() :
	protocols::moves::Mover(),
	center_(0, 0, 0),
	normal_(0, 0, 1)
{}

/// @brief Copy Constructor
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
	return ( protocols::moves::MoverOP( new MPDockingMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
MPDockingMover::fresh_instance() const {
	return protocols::moves::MoverOP( new MPDockingMover() );
}

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this Mover (MPDockingMover)
std::string
MPDockingMover::get_name() const {
	return "MPDockingMover";
}


/// @brief Add membrane components to the pose, then dock proteins along
///			the flexible jump
void MPDockingMover::apply( Pose & pose ) {
	
	using namespace core::conformation::membrane;
	
	// calling setup
	setup();
	
	// read in native pose
	read_native( pose );
	
	// assuming that protein 1 is fixed in the membrane!!!
	// add membrane VRT, call AddMembraneMover
	TR << "adding MEM" << std::endl;
	add_membrane_mover_->apply( pose );
	
	// creating foldtree from pose
	pose.fold_tree().show(std::cout);
	core::kinematics::FoldTree foldtree = pose.fold_tree();
	
	// reorder only reorders, but does not rename jump edges
	foldtree.reorder( pose.conformation().membrane_info()->membrane_rsd_num() );
	pose.fold_tree( foldtree );

	// show foldtree
	TR << "foldtree reordered" << std::endl;
	pose.fold_tree().show(std::cout);
	
	// attach Pymol observer
	TR << "test print" << std::endl;
	TR << "attach Pymol observer" << std::endl;
	//	protocols::moves::AddPyMolObserver( pose );
	
	// run docking protocol (low-res and high-res)
	TR << "calling docking protocol" << std::endl;
	docking_protocol_->apply( pose );

}

////////////////////////////////////////////////////////////////////////////////
// setup docking protocol
void MPDockingMover::setup(){

	using namespace protocols::membrane;
	using namespace protocols::docking;
	using namespace core::scoring;

	// set AddMembraneMover in protocol
	add_membrane_mover_ = AddMembraneMoverOP( new AddMembraneMover() );
	//	low_res_scorefxn_ = getScoreFunction();

	// create scorefunctions for lowres and highres
	ScoreFunctionOP lowres_scorefxn_ = ScoreFunctionFactory::create_score_function( "mpdocking_cen_14-7-23_no-penalties.wts" );
	ScoreFunctionOP highres_scorefxn_ = ScoreFunctionFactory::create_score_function( "mpdocking_fa_14-7-23_no-penalties.wts" );

	// create new docking protocol; both low-res and high-res
	docking_protocol_ = DockingProtocolOP( new DockingProtocol( 1, false, false, false, lowres_scorefxn_, highres_scorefxn_ ) );

	// get movable jump
	TR.Debug << "movable jumps: " << to_string(docking_protocol_->movable_jumps()) << std::endl;

	// set kT in docking protocol
	kT_ = 1;
}

////////////////////////////////////////////////////////////////////////////////
// read in native flags for rmsd calculation
void MPDockingMover::read_native( const Pose & pose ){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// initialize native
	core::pose::PoseOP native;
	
	// if native flag given, set native from flag, otherwise from pose
	if ( option[OptionKeys::in::file::native].user() ){
		native = core::import_pose::pose_from_pdb(option[OptionKeys::in::file::native].value_string() );
	}
	else {
		native = PoseOP( new Pose( pose ) );
	}
	
	// add membrane to native to have equal number of atoms for rmsd calculation
	add_membrane_mover_->apply( *native );
	
	// set native in docking protocol
	docking_protocol_->set_native_pose( native );
}

} // membrane
} // docking
} // protocols

#endif // INCLUDED_protocols_docking_membrane_MPDockingMover_cc
