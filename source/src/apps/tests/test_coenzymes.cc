// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    apps/tests/test_explicit_membrane_input.cc
/// @brief   Application source code for testing loading of a variety of small-molecule cofactors.
/// @author  Labonte  <JWLabonte@jhu.edu>


// Project headers
#include <devel/init.hh>

#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/minimization_packing/MinMover.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/file/FileName.hh>

// C++ headers
#include <iostream>


using namespace std;
using namespace utility;
using namespace core;
using namespace conformation;
using namespace pose;
using namespace import_pose;
using namespace scoring;


// Constants
int const SUCCESS( 0 );
int const FAILURE( -1 );

string const INPUT_PATH( "input/" );
string const OUTPUT_PATH( "output/" );


void
test_cofactor( Pose & cofactor, string const & name )
{
	cout << "---------------------------------------------------------------------------------------------" << endl;
	cout << "Importing " << name << ':' << endl;
	pose_from_file( cofactor, INPUT_PATH + name + ".pdb", PDB_file );

	cout << endl << cofactor << endl;

	file::FileName filename( name );
	filename.path( OUTPUT_PATH );
	cout << "Writing file: " << filename.ext( "pdb" ) << endl;
	cofactor.dump_pdb( filename.ext( "pdb" ) );

	cout << endl << "Residue Info:" << endl;
	Size const n_res( cofactor.size() );
	for ( core::uint i( 1 ); i <= n_res; ++i ) {
		Residue const & res( cofactor.residue( i ) );
		cout << "PDB ID: " << cofactor.pdb_info()->pose2pdb( i ) << ": ";
		res.show( cout, true );  // Show verbose output, including atomic details.
		cout << endl << endl;
	}

	// Set up ScoreFunction.
	ScoreFunctionOP sf( get_score_function() );

	cout << "Initial score: " << ( *sf )( cofactor ) << endl;

	// Set up MinMover.
	kinematics::MoveMapOP mm( new kinematics::MoveMap );
	mm->set_bb( true );
	mm->set_chi( true );
	mm->set_branches( true );
	protocols::minimization_packing::MinMover minimizer;
	minimizer.movemap( mm );
	minimizer.score_function( sf );

	minimizer.apply( cofactor );

	cout << "Score after minimization: " << ( *sf )( cofactor ) << endl;

	cofactor.dump_pdb( filename.ext( "minimized.pdb" ) );
}


int
main( int argc, char *argv[] )
{
	try {
		// Initialize core.
		devel::init( argc, argv );

		// Declare variables.
		Pose ATP, UDP_GlcNAc, SAM, SAH;
		//Pose ATP, UDP_GlcNAc, SAM, SAH, flavin1 /*, flavin2, flavin3*/ ;

		test_cofactor( ATP, "ATP" );
		test_cofactor( UDP_GlcNAc, "UDP-GlcNAc" );
		test_cofactor( SAM, "SAM" );
		test_cofactor( SAH, "SAH" );
		//test_cofactor( flavin1, "flavin_mononucleotide" );
		//test_cofactor( flavin2, "oxidized_flavin" );
		//test_cofactor( flavin3, "reduced_flavin" );

		cout << "---------------------------------------------------------------------------------------------" << endl;

	} catch (excn::Exception const & e ) {
		e.display();
		return FAILURE;
	}
	return SUCCESS;
}
