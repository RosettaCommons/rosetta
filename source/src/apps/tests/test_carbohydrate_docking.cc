// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    apps/tests/test_carbohydrate_docking.cc
/// @brief   Application source code for testing methods relevant to protein-glycoligand docking.
/// @author  Morgan L. Nance <morganlnance@gmail.com>


// Project headers
#include <devel/init.hh>

#include <core/types.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pose/util.hh>
#include <core/pose/carbohydrates/util.hh>
#include <protocols/docking/util.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <iostream>


using namespace core::import_pose;


// Constants
int const SUCCESS( 0 );
int const FAILURE( -1 );

std::string const PATH( "input/" );


void
setup_docking_ft_with_randomization( core::pose::Pose & prot_glyc_complex,
	std::string const & interface /*A_X*/ )
{
	std::cout << " FoldTree before setup_foldtree:" << std::endl;
	std::cout << prot_glyc_complex.fold_tree() << std::endl;

	// First, show what the docking FoldTree would look like using
	// the default behavior of setup_foldtree
	protocols::docking::DockJumps movable_jumps;
	protocols::docking::setup_foldtree(prot_glyc_complex,
		interface,
		movable_jumps,
		false /*randomize Jump residue of partner2*/);
	std::cout << " Docking FoldTree using default behavior (not desired):" << std::endl;
	std::cout << prot_glyc_complex.fold_tree() << std::endl;

	// Then, show the docking FoldTree using a random residue
	// from partner2 as the Jump residue
	movable_jumps.clear();
	protocols::docking::setup_foldtree(prot_glyc_complex,
		interface,
		movable_jumps,
		true /*randomize Jump residue of partner2*/);
	std::cout << " Docking FoldTree using random partner2 Jump residue (desired):" << std::endl;
	std::cout << prot_glyc_complex.fold_tree() << std::endl;

	return void();
}


int
main( int argc, char *argv[] )
{
	try {
		// Initialize core.
		devel::init( argc, argv );

		// Declare Pose variables.
		core::pose::Pose tripepA_hexaGlcX;
		std::string const & tripepA_hexaGlcX_intf = "A_X";

		// Begin testing of setup_foldtree method. /////////////////////////////////
		std::cout << "Importing tripepA_hexaGlcX.pdb" << std::endl;
		std::cout << "Using interface across chains: " <<
			tripepA_hexaGlcX_intf << std::endl;

		pose_from_file( tripepA_hexaGlcX, PATH + "tripepA_hexaGlcX.pdb", PDB_file );
		setup_docking_ft_with_randomization(tripepA_hexaGlcX, tripepA_hexaGlcX_intf);

		std::cout << std::endl << "-------------------------------------------------------------------------------" << std::endl;

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return FAILURE;
	}
	return SUCCESS;
}
