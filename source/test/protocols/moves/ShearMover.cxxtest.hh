// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/SmallMover.cxxtest.hh
/// @brief  test suite for protocols::simple_moves::SmallMover.cc
/// @author Monica Berrondo


// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <core/types.hh>

// Unit headers
#include <protocols/simple_moves/BackboneMover.hh>

// Package headers
#include <core/kinematics/MoveMap.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <core/id/AtomID_Mask.hh>
#include <utility/vector1.hh>

static basic::Tracer TR("protocols.moves.ShearMover.cxxtest");

// using declarations
using namespace core;
using namespace core::pose;
using namespace protocols::moves;

///////////////////////////////////////////////////////////////////////////
/// @name ShearMoverTest
/// @brief: test a single small move using all default values
/// @details Default values used are 0, 5, 6 (H,E,L)
///						nmoves=1, all residues movable
///  					The values for the new phi and psi torsion angles
///						were checked against values calculated by hand using
///						4 significant figures for the calculations.
///						Residue number 77 is the residue that is picked
///						The new phi angle is -74.4657 (assuming only L is being used)
///						The new psi angle is -38.8769 (for resnum-1)
///
/// @author Monica Berrondo October 09 2007
///////////////////////////////////////////////////////////////////////////
class ShearMoverTest : public CxxTest::TestSuite {

public:
	chemical::ResidueTypeSetCAP residue_set;

	PoseOP the_pose;
	ShearMoverTest() {}

	void setUp() {
		core_init();
		residue_set = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );

		the_pose = PoseOP( new Pose );
		core::import_pose::pose_from_pdb( *the_pose, "protocols/moves/test_in.pdb" );

		core::init::init_random_generators(1000, "mt19937");
	}

	void tearDown() {
		the_pose.reset();
	}

	void test_OneShearMover() {
		core::Real correct_phi ( 114.054747947984 );
		core::Real correct_psi ( 98.8718831439634 );
		// make the move
		protocols::simple_moves::ShearMover mover; // create a default small mover
		mover.apply( *the_pose );

		core::Real phi = mover.new_phi();
		core::Real psi = mover.new_psi();
		// TR.precision( 15 );
		TR << "phis and psis: " << phi << ' ' << psi << std::endl;

		// compare to the correct answer
		float const TOLERATED_ERROR = 0.0001;

		TS_ASSERT_DELTA( phi, correct_phi, TOLERATED_ERROR );
		TS_ASSERT_DELTA( psi, correct_psi, TOLERATED_ERROR );
		TR << "test_OneShearMover completed!! " << std::endl;
	}

	void test_ShearMover() {
		// this is assuming that the residue is of type L so that the max_angle is 6.0
		const int size = 12; ///< size of array to store numbers
		const int N = 1000; ///< number of trials
		const int resnum = 77; ///< residue number to change phis and psis
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
		core::Real phi = mover.new_phi();
		core::Real psi = mover.new_psi();
		// store the phi and psi as the original ones to go back to
		core::Real orig_phi = phi;
		core::Real orig_psi = psi;

		// make the move 100 times
		for ( int i=0; i<N; ++i ) {
			mover.apply( *the_pose );
			phi = mover.new_phi();
			psi = mover.new_psi();
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
		for(unsigned int i=0; i<phi_V.size(); i++) {
			if( min_c_phi > phi_V[i] ) min_c_phi = phi_V[i];
			if( max_c_phi < phi_V[i] ) max_c_phi = phi_V[i];
			if( min_c_psi > psi_V[i] ) min_c_psi = psi_V[i];
			if( max_c_psi < psi_V[i] ) max_c_psi = psi_V[i];
		}
		double phi_rate = double(min_c_phi)/max_c_phi;
		TR << min_c_phi << " ";
		TR << max_c_phi << " ";
		TR << phi_rate << " ";
		TS_ASSERT_LESS_THAN(0.4, phi_rate);
		double psi_rate = double(min_c_psi)/max_c_psi;
		TR << min_c_psi << " ";
		TR << max_c_psi << " ";
		TR << psi_rate << " ";
		TS_ASSERT_LESS_THAN(0.4, psi_rate);
		TR << "test_ShearMover completed!! " << std::endl;
	}
};

