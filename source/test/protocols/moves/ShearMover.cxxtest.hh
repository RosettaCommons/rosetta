// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/moves/SmallMover.cxxtest.hh
/// @brief  test suite for protocols::simple_moves::SmallMover.cc


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Unit headers
#include <protocols/simple_moves/BackboneMover.hh>

// Project headers
#include <core/types.hh>
#include <core/id/AtomID_Mask.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

// Utility header
#include <utility/vector1.hh>

// Basic headers
#include <basic/basic.hh>
#include <basic/Tracer.hh>


static basic::Tracer TR( "protocols.moves.ShearMover.cxxtest" );

// using declarations
using namespace core;
using namespace core::pose;
using namespace protocols::moves;


class ShearMoverTest : public CxxTest::TestSuite {

public:
	chemical::ResidueTypeSetCAP residue_set;

	PoseOP the_pose;
	PoseOP the_sugar_pose;
	ShearMoverTest() {}

	void setUp() {
		core_init_with_additional_options( "-include_sugars" );
		residue_set = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );

		the_pose = PoseOP( new Pose );
		core::import_pose::pose_from_file( *the_pose, "protocols/moves/test_in.pdb" , core::import_pose::PDB_file);

		the_sugar_pose = PoseOP( new Pose );
		make_pose_from_saccharide_sequence( *the_sugar_pose,
			"b-D-Glcp-(1->4)-b-D-Glcp-(1->4)-[b-D-Glcp-(1->4)-[b-D-Glcp-(1->6)]-b-D-Glcp-(1->6)]-"
			"b-D-Glcp-(1->4)-b-D-Glcp-(1->4)-b-D-Glcp" );

		core::init::init_random_generators(1000, "mt19937");
	}

	void tearDown() {
		the_pose.reset();
	}

	/// @author Labonte <JWLabonte@jhu.edu>
	void test_one_shear_move() {
		protocols::simple_moves::ShearMover mover;

		TR << "Testing that the ShearMover behaves properly with a simple peptide." << std::endl;

		// Force the 76th and 77th residues of this pose to move.
		core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
		mm->set_bb( false );
		mm->set_bb( 76, true );
		mm->set_bb( 77, true );
		mover.movemap( mm );

		Angle const initial_phi( the_pose->phi( 77 ) ), initial_psi( the_pose->psi( 76 ) );
		mover.apply( *the_pose );

		TS_ASSERT( the_pose->phi( 77 ) != initial_phi );
		TS_ASSERT( the_pose->psi( 76 ) != initial_psi );
		TS_ASSERT_DELTA( basic::periodic_range( the_pose->phi( 77 ) - initial_phi, 360 ),
			basic::periodic_range( initial_psi - the_pose->psi( 76 ), 360 ), 0.00001 );
	}

	/// @author Monica Berrondo
	void test_ShearMover() {
		// this is assuming that the residue is of type L so that the max_angle is 6.0
		const int size = 12; //< size of array to store numbers
		const int N = 1000; //< number of trials
		const int resnum = 77; //< residue number to change phis and psis
		// create array and fill it
		std::vector<int> phi_V(size), psi_V(size);
		for ( unsigned int i=0; i<phi_V.size(); ++i )  {
			phi_V[i] = 0;
			psi_V[i] = 0;
		}

		// set up the mover to only change on residue and loop 100 times
		protocols::simple_moves::ShearMover mover; // create a default small mover
		core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
		movemap->set_bb( false );
		movemap->set_bb( resnum, true );
		movemap->set_bb( resnum-1, true );
		mover.movemap( movemap );

		// make an initial move and get the phi and psis
		mover.apply( *the_pose );
		core::Real phi = the_pose->phi( resnum );
		core::Real psi = the_pose->psi( resnum - 1 );
		// store the phi and psi as the original ones to go back to
		core::Real orig_phi = phi;
		core::Real orig_psi = psi;

		// make the move 100 times
		for ( int i=0; i<N; ++i ) {
			mover.apply( *the_pose );
			phi = the_pose->phi( resnum );
			psi = the_pose->psi( resnum - 1 );
			int phi_ind = int( phi - orig_phi + int(size/2) );
			int psi_ind = int( psi - orig_psi + int(size/2) );
			phi_V[phi_ind] += 1;
			psi_V[psi_ind] += 1;

			the_pose->set_phi( resnum, orig_phi );
			the_pose->set_psi( resnum-1, orig_psi );
		}

		for ( unsigned int i=0; i<phi_V.size(); ++i )  {
			TR << "phi of " << i << " : " << phi_V[i];
			TR << " / psi of " << i << " : " << psi_V[i] << std::endl;
		}

		int min_c_phi = N;
		int max_c_phi = 0;
		int min_c_psi = N;
		int max_c_psi = 0;
		for ( unsigned int i=0; i<phi_V.size(); i++ ) {
			if ( min_c_phi > phi_V[i] ) min_c_phi = phi_V[i];
			if ( max_c_phi < phi_V[i] ) max_c_phi = phi_V[i];
			if ( min_c_psi > psi_V[i] ) min_c_psi = psi_V[i];
			if ( max_c_psi < psi_V[i] ) max_c_psi = psi_V[i];
		}
		double phi_rate = double(min_c_phi)/max_c_phi;
		TR << min_c_phi << " ";
		TR << max_c_phi << " ";
		TR << phi_rate << " ";
		TS_ASSERT_LESS_THAN(/*0.4*/0.35, phi_rate); // decreases with beta_nov15
		double psi_rate = double(min_c_psi)/max_c_psi;
		TR << min_c_psi << " ";
		TR << max_c_psi << " ";
		TR << psi_rate << " ";
		TS_ASSERT_LESS_THAN(/*0.4*/0.35, psi_rate); // decreases with beta_nov15
		TR << "test_ShearMover completed!! " << std::endl;
	}

	/// @author Labonte <JWLabonte@jhu.edu>
	void test_ShearMover_w_sugar() {
		protocols::simple_moves::ShearMover mover;
		core::select::residue_selector::ResidueIndexSelectorOP selector(
			new core::select::residue_selector::ResidueIndexSelector );
		core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
		mm->set_bb( true );
		mm->set_chi( true );
		mm->set_branches( true );
		mover.movemap( mm );

		TR << "Testing that the ShearMover behaves properly with a complicated sugar." << std::endl;

		// Make a shear move at residue 2, an internal residue with two glycosidic torsions.
		// Since residue 1 is a lower terminus and has no glycosidic torsions, this should force counter moves from
		// residue 3.
		selector->set_index( "2" );
		mover.selector( selector );

		Angle initial_phi2( the_sugar_pose->phi( 2 ) ), initial_psi2( the_sugar_pose->psi( 2 ) ),
			initial_phi3( the_sugar_pose->phi( 3 ) ), initial_psi3( the_sugar_pose->psi( 3 ) );
		mover.apply( *the_sugar_pose );

		if ( the_sugar_pose->phi( 2 ) != initial_phi2 ) {
			// If phi2 was randomly selected, psi3 should be the counter move, since the 4C1 ring of all-equatorial
			// glucose holds them nearly parallel.
			TS_ASSERT( the_sugar_pose->psi( 2 ) == initial_psi2 );
			TS_ASSERT( the_sugar_pose->phi( 3 ) == initial_phi3 );
			TS_ASSERT( the_sugar_pose->psi( 3 ) != initial_psi3 );
			TS_ASSERT_DELTA( basic::periodic_range( the_sugar_pose->phi( 2 ) - initial_phi2, 360 ),
				basic::periodic_range( initial_psi3 - the_sugar_pose->psi( 3 ), 360 ), 0.00001 );
		} else if ( the_sugar_pose->psi( 2 ) != initial_psi2 ) {
			// If psi2 was randomly selected, we can't know for certain what the best counter move is without more
			// information.
			TS_ASSERT( true );
		} else {
			// If neither phi2 nor psi2 moved, the Mover failed to work properly.
			TS_ASSERT( false );
		}


		// Make a shear move at residue 5, an upper terminus residue with two glycosidic torsions.
		// This should force counter moves from residue 4.
		selector->set_index( "5" );
		mover.selector( selector );

		Angle initial_phi5( the_sugar_pose->phi( 5 ) ), initial_psi5( the_sugar_pose->psi( 5 ) ),
			initial_phi4( the_sugar_pose->phi( 4 ) ), initial_psi4( the_sugar_pose->psi( 4 ) );
		mover.apply( *the_sugar_pose );

		if ( the_sugar_pose->phi( 5 ) != initial_phi5 ) {
			// If phi5 was randomly selected, we can't know for certain what the best counter move is without more
			// information.
			TS_ASSERT( true );
		} else if ( the_sugar_pose->psi( 5 ) != initial_psi5 ) {
			// If psi5 was randomly selected, phi4 should be the counter move, since the 4C1 ring of all-equatorial
			// glucose holds them nearly parallel.
			TS_ASSERT( the_sugar_pose->phi( 5 ) == initial_phi5 );
			TS_ASSERT( the_sugar_pose->phi( 4 ) != initial_phi4 );
			TS_ASSERT( the_sugar_pose->psi( 4 ) == initial_psi4 );
			TS_ASSERT_DELTA( basic::periodic_range( the_sugar_pose->psi( 5 ) - initial_psi5, 360 ),
				basic::periodic_range( initial_phi4 - the_sugar_pose->phi( 4 ), 360 ), 0.00001 );
		} else {
			// If neither phi5 nor psi5 moved, the Mover failed to work properly.
			TS_ASSERT( false );
		}


		// Make a shear move at residue 4, an internal residue with two glycosidic torsions.
		// Since residue 3 is a branch point, this should force counter moves from residue 5.
		selector->set_index( "4" );
		mover.selector( selector );

		initial_phi5 = the_sugar_pose->phi( 5 );
		initial_psi5 = the_sugar_pose->psi( 5 );
		initial_phi4 = the_sugar_pose->phi( 4 );
		initial_psi4 = the_sugar_pose->psi( 4 );
		mover.apply( *the_sugar_pose );

		if ( the_sugar_pose->phi( 4 ) != initial_phi4 ) {
			// If phi4 was randomly selected, psi5 should be the counter move, since the 4C1 ring of all-equatorial
			// glucose holds them nearly parallel.
			TS_ASSERT( the_sugar_pose->psi( 4 ) == initial_psi4 );
			TS_ASSERT( the_sugar_pose->phi( 5 ) == initial_phi5 );
			TS_ASSERT( the_sugar_pose->psi( 5 ) != initial_psi5 );
			TS_ASSERT_DELTA( basic::periodic_range( the_sugar_pose->phi( 4 ) - initial_phi4, 360 ),
				basic::periodic_range( initial_psi5 - the_sugar_pose->psi( 5 ), 360 ), 0.00001 );
		} else if ( the_sugar_pose->psi( 4 ) != initial_psi4 ) {
			// If psi4 was randomly selected, we can't know for certain what the best counter move is without more
			// information.
			TS_ASSERT( true );
		} else {
			// If neither phi4 nor psi4 moved, the Mover failed to work properly.
			TS_ASSERT( false );
		}


		// Make a shear move at residue 3, a branch point with two glycosidic torsions.
		// This should force counter moves from residue 2.
		selector->set_index( "3" );
		mover.selector( selector );

		initial_phi3 = the_sugar_pose->phi( 3 );
		initial_psi3 = the_sugar_pose->psi( 3 );
		initial_phi2 = the_sugar_pose->phi( 2 );
		initial_psi2 = the_sugar_pose->psi( 2 );
		mover.apply( *the_sugar_pose );

		if ( the_sugar_pose->phi( 3 ) != initial_phi3 ) {
			// If phi3 was randomly selected, we can't know for certain what the best counter move is without more
			// information.
			TS_ASSERT( true );
		} else if ( the_sugar_pose->psi( 3 ) != initial_psi3 ) {
			// If psi3 was randomly selected, phi2 should be the counter move, since the 4C1 ring of all-equatorial
			// glucose holds them nearly parallel.
			TS_ASSERT( the_sugar_pose->phi( 3 ) == initial_phi3 );
			TS_ASSERT( the_sugar_pose->phi( 2 ) != initial_phi2 );
			TS_ASSERT( the_sugar_pose->psi( 2 ) == initial_psi2 );
			TS_ASSERT_DELTA( basic::periodic_range( the_sugar_pose->psi( 3 ) - initial_psi3, 360 ),
				basic::periodic_range( initial_phi2 - the_sugar_pose->phi( 2 ), 360 ), 0.00001 );
		} else {
			// If neither phi3 nor psi3 moved, the Mover failed to work properly.
			TS_ASSERT( false );
		}


		// Make a shear move at residue 6, an internal residue with three glycosidic torsions.
		// Since both residue 5 and residue 6 are branch points, only torsions on residue 6 can move, and since there
		// are three glycosidic torsions, only phi6 and omega6 can counter each other.
		selector->set_index( "6" );
		mover.selector( selector );

		Angle initial_phi6( the_sugar_pose->phi( 6 ) ), initial_psi6( the_sugar_pose->psi( 6 ) ),
			initial_omega6( the_sugar_pose->omega( 6 ) );
		mover.apply( *the_sugar_pose );

		if ( the_sugar_pose->phi( 6 ) != initial_phi6 ) {
			// If phi6 was randomly selected, omega6 should be the counter move, since its the only option.  (psi6 is
			// adjacent and could not possibly be near-parallel.
			TS_ASSERT( the_sugar_pose->psi( 6 ) == initial_psi6 );
			TS_ASSERT( the_sugar_pose->omega( 6 ) != initial_omega6 );
			TS_ASSERT_DELTA( basic::periodic_range( the_sugar_pose->phi( 6 ) - initial_phi6, 360 ),
				basic::periodic_range( initial_omega6 - the_sugar_pose->omega( 6 ), 360 ), 0.00001 );
		} else if ( the_sugar_pose->omega( 6 ) != initial_omega6 ) {
			// If omega6 was randomly selected, phi6 should be the counter move, since its the only option.  (psi6 is
			// adjacent and could not possibly be near-parallel.
			TS_ASSERT( the_sugar_pose->phi( 6 ) != initial_phi6 );
			TS_ASSERT( the_sugar_pose->psi( 6 ) == initial_psi6 );
			TS_ASSERT_DELTA( basic::periodic_range( the_sugar_pose->omega( 6 ) - initial_omega6, 360 ),
				basic::periodic_range( initial_phi6 - the_sugar_pose->phi( 6 ), 360 ), 0.00001 );
		} else {
			// If neither phi6 nor omega6 moved, the Mover failed to work properly.
			TS_ASSERT( false );
		}
	}
};
