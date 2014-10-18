// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		apps/pilot/membrane/membrane_symdocking.cc
///
/// @brief		Membrane Framework Application: Symmetric Protein-Protein Docking in Membranes
/// @details	Below is a first-pass at plugging in the membrane framework with symmetry.
///				Uses membrane representation, energ functions, and current tools for Rosetta
///				design for computing ddGs.
///				Last Modified: 8/31/14
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <protocols/moves/Mover.hh> 

#include <protocols/symmetric_docking/SymDockProtocol.hh> 
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>

#include <protocols/membrane/symmetry/SetupForMembraneSymmetry.hh> 

#include <protocols/membrane/AddMembraneMover.hh> 

// Project Headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>

#include <core/scoring/ScoreFunction.hh> 
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/conformation/Conformation.hh> 
#include <core/conformation/membrane/MembraneInfo.hh>

// symm res - start
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/kinematics/FoldTree.hh>

// sym res - end

#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/pointer/owning_ptr.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <iostream>

using namespace core::pose;
using namespace core::scoring;
using namespace protocols::membrane;
using namespace protocols::simple_moves::symmetry;
using namespace protocols::symmetric_docking;

static basic::Tracer TR( "apps.pilot.ralford.mp_symdocking" );

/// @brief Mover for Symmetric Docking in Membranes
class MPSymDockMover : public protocols::moves::Mover {

public:
	
	/// @brief Default Constructor
	MPSymDockMover() : Mover() {}
	
	/// @brief Required get_name method
	std::string get_name() const { return "MPSymDockMover"; }
	
	/// @brief Symmetric Docking in the Membrane
	/// @details Setup pose for symmetry, add the membrane components to the total pose,
	/// and then perform symmetric docking
	void apply( Pose & pose ) {
	
		using namespace protocols::membrane::symmetry; 

		// Setting up the membrane energy functions (low & high resolution)
		ScoreFunctionOP sfxn_low = ScoreFunctionFactory::create_score_function( "mpsymdocking_lowres" );
		ScoreFunctionOP sfxn_high = ScoreFunctionFactory::create_score_function( "mpsymdocking_hires" );
	

		//TR << "Creating setup symm" << std::endl;

		// Setup the pose for Symmetry
		// If pose is a membrane enabled pose, choose the special membrane residue as the root
		//SetupForSymmetryMoverOP setup_symm( new SetupForSymmetryMover() );
		//setup_symm->apply( pose );

		//TR << "Creating add memb" << std::endl;
		
		// Add membrane to pose
		//AddMembraneMoverOP add_memb( new AddMembraneMover() );
		//add_memb->apply( pose );


		SetupForMembraneSymmetryOP setup_memb( new SetupForMembraneSymmetry() ); 
		setup_memb->apply( pose ); 

		// Rematching symmetry effects
		core::Vector center(0, 0, 0); 
		core::Vector normal(0, 0, 1); 
		pose.conformation().update_membrane_position(center, normal);

		TR << pose.conformation().membrane_info()->membrane_rsd_num() << std::endl;

		
		//// Everything in this box can go into a setup for membrane symmetry mover ////////
		
		// Setup Membranes for Symmetry
		// ---------------------------------------------------------------------
		// 1. Setup new conformaiton pointer
		// 2. Add master membrane subunit to the root of the symmetric fold tree
		// 3. Inform MembraneInfo that the master membrane residue for scoring has changed, as
		// well as the jump of kinematic control
		// 4. compute summed membrane from child residues and set this new membrane
		// as the master
		// (5. some method - compute symmetric embedding in the membrane?)
		
		//// Closing the box /////////////////////////////////////////////////////////////////
		
		// Go do some docking!
		SymDockProtocolOP symdock( new SymDockProtocol( true, false, false, sfxn_low, sfxn_high ) );
		symdock->apply( pose );
		
	}
};

// typedefs
typedef utility::pointer::shared_ptr< MPSymDockMover > MPSymDockMoverOP;

/// @brief Main method
int
main( int argc, char * argv [] )
{
	using namespace protocols::jd2;
	
	try {
		
		// Devel init factories
		devel::init(argc, argv);
		
		// Register JD2 options
		protocols::jd2::register_options();
		
		// Setup Membrane Symdocking & go!
		MPSymDockMoverOP symdock( new MPSymDockMover() ); 
		protocols::jd2::JobDistributor::get_instance()->go( symdock );

		return 0;
		
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
