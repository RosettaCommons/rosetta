// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file 		src/apps/pilot/ralford/membrane_relax2.cc
///
/// @brief 		Membrane Protein Structure - FastRelax Protocol (Structure Refinement)
/// @details	Custom version of FastRelax - sample local conformational change of protein
///	     	    structure in the membrane. Steps: (1) membrane insertion and (2) execting cycles of repulsive 
///  			upweighting, repack and minimization of structure with respect to the membrane.
///
/// 			Last Modified: 6/28/14
///				Note: Part of the Membrane Framework Applications Suite
///				Versioning: 1.0
///
/// @author 	Rebecca Alford (rfalford12@gmail.com)

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <protocols/membrane/AddMembraneMover.hh> 
#include <protocols/membrane/MembranePositionFromTopologyMover.hh> 

// Project Headers
#include <protocols/jd2/JobDistributor.hh> 
#include <protocols/jd2/util.hh>

#include <core/pose/Pose.hh> 
#include <core/conformation/Conformation.hh>

#include <core/scoring/ScoreFunction.hh> 
#include <core/scoring/ScoreFunctionFactory.hh>  

#include <core/kinematics/MoveMap.hh> 

#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/AtomTreeMinimizer.hh> 

#include <protocols/simple_moves/MinMover.hh> 

// Utility Headers
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

#include <numeric/random/random.hh> 
#include <utility/pointer/owning_ptr.hh> 

#include <numeric/conversions.hh> 

#include <basic/Tracer.hh> 

// C++ headers
#include <iostream>

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.ralford.membrane_relax2" );

using namespace core; 
using namespace protocols::moves;

/// @brief Membrane Transform Mover: Transform a pose into a new mover
class MembraneRelaxMover : public Mover {

public: 

	MembraneRelaxMover() {}

	/// @brief Get Mover Name
	std::string get_name() const { return "MembraneRelaxMover"; }

	/// @brief Apply Membrane Relax
	void apply( Pose & pose ) {

		using namespace core::scoring; 
		using namespace core::kinematics; 
		using namespace core::optimization; 
		using namespace protocols::membrane;

		// Add a Membrane
		AddMembraneMoverOP add_memb = new AddMembraneMover(); 
		add_memb->apply(pose); 

		// Initialize Membrane Position
		MembranePositionFromTopologyMoverOP init_memb = new MembranePositionFromTopologyMover(); 
		init_memb->apply( pose ); 

		// Setup Membrane Energy Function as is
		ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( "mpframework_fa_2007" );

		// Setup a MoveMap
		MoveMapOP movemap = new MoveMap();
		movemap->set_bb( 1, 222 );
		movemap->set_chi( true ); 
		movemap->set_jump( true );

		// Setup minimizer options
		core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_nonmonotone", 0.1, true );
		min_opts.max_iter( 200 ); 

		// Setup Minimizer
		AtomTreeMinimizer at_min;
		(*sfxn)(pose);

		// Grab jump CA anchors
		core::Vector atom1a = pose.residue( 1 ).xyz( 2 ); 
		core::Vector atom2a = pose.residue( 223 ).xyz( 2 ); 

		TR << "Initial Membrnae Position: " << std::endl;
		show_membrane_position( pose ); 

		TR << "Initial Distance: " << atom1a.distance( atom2a ) << std::endl;

		// Minimize in Rounds
		//cart_min.run( pose, *movemap, *sfxn, min_opts ); 
		at_min.run( pose, *movemap, *sfxn, min_opts ); 

		TR << "Final Membrane Position: " << std::endl;
		show_membrane_position( pose ); 

		// Compute Jump Distance
		core::Vector atom1b = pose.residue( 1 ).xyz( 2 ); 
		core::Vector atom2b = pose.residue( 223 ).xyz( 2 ); 

		TR << "Final distance: " << atom1b.distance( atom2b ) << std::endl;
	}

	void
	show_membrane_position( Pose & pose ) {

		Vector center( pose.conformation().membrane_info()->membrane_center(pose.conformation()) ); 
		Vector normal( pose.conformation().membrane_info()->membrane_normal(pose.conformation()) );

		TR << "Membrane Center: " << center.x() << " " << center.y() << " " << center.z() << std::endl;
		TR << "Membrane Normal: " << normal.x() << " " << normal.y() << " " << normal.z() << std::endl;
	}
};

typedef utility::pointer::owning_ptr< MembraneRelaxMover > MembraneRelaxMoverOP; 
typedef utility::pointer::owning_ptr< MembraneRelaxMover const > MembraneRelaxMoverCOP; 

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
		MembraneRelaxMoverOP mprelax = new MembraneRelaxMover();
		protocols::jd2::JobDistributor::get_instance()->go( mprelax );

		return 0; 

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

