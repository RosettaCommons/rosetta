// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    apps/tests/test_explicit_membrane_input.cc
/// @brief   Application source code for testing loading of lipids within the context of a small surface.
/// @author  Labonte  <JWLabonte@jhu.edu>


// Project headers
#include <devel/init.hh>

#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/simple_moves/MinMover.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <iostream>


using namespace std;
using namespace core;
using namespace pose;
using namespace import_pose;
using namespace conformation;
using namespace scoring;


string const INPUT_PATH( "input/" );
string const OUTPUT_PATH( "output/" );


int
main( int argc, char *argv[] )
{
	try {
		// Initialize core.
		devel::init( argc, argv );

		// Declare variables.
		Pose lipid_surface;


		cout << "Importing lipid monolayer with 99 lipids:" << endl;

		pose_from_file( lipid_surface, INPUT_PATH + "lipid_surface.pdb", PDB_file );

		cout << endl << lipid_surface << endl;
 
		cout << endl << "Residue Info:" << endl;
		Size const n_res( 8 );  // Only output the first two lipids, each composed of 4 residues.
		for ( core::uint i( 1 ); i <= n_res; ++i ) {
			Residue const & res( lipid_surface.residue( i ) );
			cout << "PDB ID: " << lipid_surface.pdb_info()->pose2pdb( i ) << ": ";
			res.show( cout, true );  // Show verbose output, including atomic details.
			cout << endl << endl;
		}

		lipid_surface.dump_pdb( OUTPUT_PATH + "lipid_surface_initial.pdb" );

		// Set up ScoreFunction.
		ScoreFunctionOP sf( get_score_function() );

		cout << "Initial score of membrane complex: " << ( *sf )( lipid_surface ) << endl;

		// Set up MinMover.
		kinematics::MoveMapOP mm( new kinematics::MoveMap );
		mm->set_bb( true );
		mm->set_chi( true );
		mm->set_branches( true );
		protocols::simple_moves::MinMover minimizer;
		minimizer.movemap( mm );
		minimizer.score_function( sf );

		minimizer.apply( lipid_surface );

		cout << "Score after minimization: " << ( *sf )( lipid_surface ) << endl;

		lipid_surface.dump_pdb( OUTPUT_PATH + "lipid_surface_minimized.pdb" );

		cout << "---------------------------------------------------------------------------------------------" << endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		cerr << "Caught exception: " << e.msg() << endl;
		return -1;
	}
	return 0;
}
