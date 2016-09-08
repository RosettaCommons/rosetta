// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file 		src/apps/pilot/ralford/membrane_sampling.cc
///
/// @brief 		Pilot Application for Testing Membrane & Non membrane movers
///
/// @author 	Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified: 7/10/14

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <protocols/moves/Mover.hh>
#include <protocols/rigid/RigidBodyMover.hh> 

#include <protocols/jd2/JobDistributor.hh> 
#include <protocols/jd2/util.hh>

#include <protocols/membrane/AddMembraneMover.hh> 
#include <protocols/membrane/MembranePositionFromTopologyMover.hh> 

#include <core/scoring/ScoreFunction.hh> 
#include <core/scoring/ScoreFunctionFactory.hh> 

#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh> 
#include <core/optimization/MinimizerOptions.hh> 

#include <core/conformation/Conformation.hh> 
#include <core/conformation/membrane/MembraneInfo.hh> 

#include <core/pose/Pose.hh> 
#include <core/types.hh> 

// Utility Headers
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

#include <numeric/xyz.functions.hh> 
#include <numeric/xyzVector.hh> 

#include <numeric/random/random.hh> 

#include <utility/pointer/owning_ptr.hh> 

#include <basic/Tracer.hh> 

using namespace protocols::moves;

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.ralford.membrane_sampling" );

/// @brief View membrane planes based on normal/center
class MembraneSampleMover : public Mover {

public: 

	MembraneSampleMover() {}

	/// @brief Get Mover Name
	std::string get_name() const { return "MembraneSampleMover"; }

	/// @brief Apply Membrane Relax
	void apply( Pose & pose ) {

		using namespace numeric; 
		using namespace core; 
		using namespace core::scoring; 
		using namespace core::kinematics;
		using namespace protocols::membrane;

		// Add a Membrane
		AddMembraneMoverOP add_memb = new AddMembraneMover(); 
		add_memb->apply( pose );

	
		// Setup the Membrane Energy Function
		ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( "mpframework_fa_2007" );

		// Set up a movemap - bb, chi and jump are all moveable
		MoveMapOP movemap = new MoveMap(); 
		movemap->set_bb( true ); 
		movemap->set_chi( true ); 
		movemap->set_jump( true ); // will need to play with this one

		// Apply the Min Mover
		core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_nonmonotone", 0.1, true );
	 	core::optimization::AtomTreeMinimizer atm;
	 	atm.run( pose, *movemap, *sfxn, min_opts ); 

		// Make 5 random phi/psi moves
	/*	for ( Size i = 1; i <= 20; ++i ) {

			core::Size randres = numeric::random::random_range( 1, pose.size() );

   		 	core::Real phi = pose.phi( randres ) + numeric::random::gaussian() * 4;
    	 	core::Real psi = pose.psi( randres ) + numeric::random::gaussian() * 4;

		    pose.set_phi( randres, phi );
    		pose.set_psi( randres, psi );

    		
		}

	*/
	}

};

typedef utility::pointer::owning_ptr< MembraneSampleMover > MembraneSampleMoverOP;
typedef utility::pointer::owning_ptr< MembraneSampleMover const > MembraneSampleMoverCOP; 

/// @brief Main method
int
main( int argc, char * argv [] )
{
	try {

		// Devel init factories
		devel::init(argc, argv);

		// Register JD2 options
		protocols::jd2::register_options();

		// Minimize with mp
		MembraneSampleMoverOP mp = new MembraneSampleMover();
		protocols::jd2::JobDistributor::get_instance()->go( mp );

		return 0; 

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
