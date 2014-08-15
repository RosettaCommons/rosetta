// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license.
// (c) The Rosetta software is developed by the contributing members of the
// (c) Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org.
// (c) Questions about this can be addressed to University of Washington UW
// (c) TechTransfer, email: license@u.washington.edu.

/// @file    test_CCD_loop_closure.cc
/// @brief   Integration test application for testing CCD loop closure.
/// @details This application...
/// ...creates polyalanine chains of varying lengths,
/// ...randomly opens the chains, and
/// ...uses CCD to close the broken loops.
/// @author  Labonte <JWLabonte@jhu.edu>

// Project headers
#include <devel/init.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/TorsionID.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh>

// Numeric header
#include <numeric/random/random.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>

// C++ header
#include <iostream>


// Using which namespaces?
using namespace std;
using namespace utility;
using namespace core;
using namespace core::id;
using namespace protocols::loops;
using namespace protocols::loops::loop_closure::ccd;
using namespace pose;
using namespace kinematics;
using namespace core::conformation;


// Construct random-number generator.
static numeric::random::RandomGenerator RG( 21 );  // the 6th triangular number


// Constants
string const PATH = "output/";


// Functions
// Get a linear, n-mer polyalanine.
PoseOP
get_n_mer_polyalanine( Size n )
{
	string sequence;
	for ( core::uint i = 1; i <= n; ++i ) {
		sequence += "A";
	}

	PoseOP pose( new Pose );
	make_pose_from_sequence( *pose, sequence, "fa_standard" );

	for ( core::uint i = 1; i < n; ++i ) {
		pose->set_omega( i, 180.0 );  // Make linear.
	}

	return pose;
}

// Get a MoveMap where only phi & psi can move, up to residue n.
MoveMapOP
get_phi_psi_mm( Size n )
{
	MoveMapOP mm( new MoveMap );
	mm->set_bb( true );
	for ( core::uint i = 1; i < n; ++i ) {
		mm->set( TorsionID( i, BB, omega_torsion ), false );
	}

	return mm;
}

// Get an appropriate loop for an n-mer, with the cut in the center and the first and last residues as anchors.
Loop
get_loop_for_n_mer( Size n )
{
	return Loop( 2, n -1 , n / 2 );
}


// Randomly open this pose's loop.
void
randomly_open_pose_loop( Pose & pose, Loop const & loop )
{
	for ( core::uint i = loop.start(); i <= loop.stop(); ++i ) {
		pose.set_phi( i, RG.uniform() * 360 );
		pose.set_psi( i, RG.uniform() * 360 );
	}
}


int
main( int argc, char *argv[] )
{
    try {
		// Initialize core.
		devel::init( argc, argv );

		// Construct the Mover.
		CCDLoopClosureMover mover;
		mover.check_rama_scores( false );
		mover.max_per_move_torsion_delta_per_residue( 180., 180., 180. );
		mover.max_total_torsion_delta_per_residue( 180., 180., 180. );
		mover.max_cycles( 500 );

		// Test closure of polyalanines of varying lengths.
		for ( Size n = 6; n <= 14; n += 2 ) {
			PoseOP pose( get_n_mer_polyalanine( n ) );
			MoveMapOP mm( get_phi_psi_mm( n ) );
			Loop const loop( get_loop_for_n_mer( n ) );

			set_single_loop_fold_tree( *pose, loop );
			add_single_cutpoint_variant( *pose, loop );
			cout << "------------------------------------------------------------------------------" << endl;
			cout << n << "-mer polyalanine: " << *pose << endl;

			cout << " Randomly opening loop..." << endl;
			randomly_open_pose_loop( *pose, loop );

			string const filename_start( PATH + "CCD_test_" + to_string( n ) + "mer_start.pdb" );
			cout << " Outputting starting structure to " << filename_start << "..." << endl;
			pose->dump_pdb( filename_start );

			cout << " Closing the loop with CCD..." << endl;
			mover.movemap( mm );
			mover.loop( loop );
			mover.apply( *pose );

			string const filename_end( PATH + "CCD_test_" + to_string( n ) + "mer_end.pdb" );
			cout << " Outputting ending structure to " << filename_end << "..." << endl;
			pose->dump_pdb( filename_end );

			cout << " Finished testing " << n << "-mer polyalanine." << endl << endl;
		}

		cout << "------------------------------------------------------------------------------" << endl;
		cout << "Done!" << endl;
	} catch ( utility::excn::EXCN_Base const & e ) {
		cerr << "Caught exception: " << e.msg() << endl;
		return -1;
	}
	return 0;
}
