// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// Copyright in the Rosetta software belongs to the developers and their institutions.
// For more information, see www.rosettacommons.org.
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/fragment/IndependentBBTorsionSRFD.cxxtest.hh
/// @brief  unit tests for IndependentBBTorsionSRFD
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/fragment/IndependentBBTorsionSRFD.hh>
#include <core/id/TorsionID.hh>
// AUTO-REMOVED #include <core/kinematics/Edge.hh>
// AUTO-REMOVED #include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>

#include <numeric/angle.functions.hh>

#include <cmath>
#include <string>

//Auto Headers
#include <core/id/AtomID_Mask.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/pose/annotated_sequence.hh>
#include <utility/fix_boinc_read.hh>
#include <utility/vector1.hh>


class IndependentBBTorsionSRFDTests : public CxxTest::TestSuite
{


public: // setup


	typedef std::string String;
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::Vector Vector;
	typedef core::pose::Pose Pose;
	typedef core::fragment::IndependentBBTorsionSRFD IndependentBBTorsionSRFD;
	typedef core::kinematics::Edge Edge;
	typedef core::kinematics::FoldTree FoldTree;
	typedef core::kinematics::MoveMap MoveMap;


	IndependentBBTorsionSRFDTests() {};


	// Shared initialization.
	void setUp() {
		core_init();
	}


	// Shared finalization.
	void tearDown() {
	}


public: // re-used methods


	/// @brief return a Pose with a continuous topology, helical geometry
	Pose helix10_pose() {
		Pose pose;
		core::pose::make_pose_from_sequence(
			pose,
			"K[LYS:NtermProteinFull]EEEEEEEEK[LYS:CtermProteinFull]",
			*core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID )
		);

		for ( core::Size i = 1, ie = pose.n_residue(); i <= ie; ++i ) {
			pose.set_secstruct( i, 'H' );
		}

		for ( Size i = 1, ie = pose.n_residue(); i <= ie; ++i ) {
			pose.set_phi( i, -60.0 );
			pose.set_psi( i, -45.0 );
			pose.set_omega( i, 180.0 );
		}

		return pose;
	}


public: // tests


	/// @brief test is_applicable() and apply()
	void test_apply() {
		using core::id::TorsionID;

		using numeric::nonnegative_principal_angle_degrees;

		Pose h10 = helix10_pose();

		IndependentBBTorsionSRFD srfd( 3, 'E', 'A' );
		srfd.set_torsion( 1, 30.0 );
		srfd.set_torsion( 2, 40.0 );
		srfd.set_torsion( 3, 50.0 );

		MoveMap movemap;

		// test empty movemap, should not be applicable
		TS_ASSERT( !srfd.is_applicable( movemap, 6 ) );
		TS_ASSERT( !srfd.apply( movemap, h10, 6 ) ); // pose should not change here except for secstruct
		TS_ASSERT_EQUALS( h10.secstruct( 6 ), 'E' ); // non-intuitive behavior on part of SecstructSRFD
		TS_ASSERT_DELTA( h10.phi( 6 ), -60.0, 0.001 );
		TS_ASSERT_DELTA( h10.psi( 6 ), -45.0, 0.001 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( h10.omega( 6 ) ), 180.0, 0.001 );

		// test locked down position at movemap bb level, should not be applicable
		movemap.set_bb( 5, false );
		TS_ASSERT( !srfd.is_applicable( movemap, 5 ) );
		TS_ASSERT( !srfd.apply( movemap, h10, 5 ) ); // pose should not change here except for secstruct
		TS_ASSERT_EQUALS( h10.secstruct( 5 ), 'E' ); // non-intuitive behavior on part of SecstructSRFD
		TS_ASSERT_DELTA( h10.phi( 5 ), -60.0, 0.001 );
		TS_ASSERT_DELTA( h10.psi( 5 ), -45.0, 0.001 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( h10.omega( 5 ) ), 180.0, 0.001 );

		// test locked down position at torsion level, should not be applicable
		movemap.clear();
		movemap.set( TorsionID( 4, core::id::BB, core::id::phi_torsion ), false );
		movemap.set( TorsionID( 4, core::id::BB, core::id::psi_torsion ), false );
		movemap.set( TorsionID( 4, core::id::BB, core::id::omega_torsion ), false );
		TS_ASSERT( !srfd.is_applicable( movemap, 4 ) );
		TS_ASSERT( !srfd.apply( movemap, h10, 4 ) ); // pose should not change here except for secstruct
		TS_ASSERT_EQUALS( h10.secstruct( 4 ), 'E' ); // non-intuitive behavior on part of SecstructSRFD
		TS_ASSERT_DELTA( h10.phi( 4 ), -60.0, 0.001 );
		TS_ASSERT_DELTA( h10.psi( 4 ), -45.0, 0.001 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( h10.omega( 4 ) ), 180.0, 0.001 );

		// open up phi
		movemap.set( TorsionID( 4, core::id::BB, core::id::phi_torsion ), true );
		TS_ASSERT( srfd.is_applicable( movemap, 4 ) );
		TS_ASSERT( srfd.apply( movemap, h10, 4 ) ); // phi changes here
		TS_ASSERT_EQUALS( h10.secstruct( 4 ), 'E' ); // non-intuitive behavior on part of SecstructSRFD
		TS_ASSERT_DELTA( h10.phi( 4 ), 30.0, 0.001 );
		TS_ASSERT_DELTA( h10.psi( 4 ), -45.0, 0.001 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( h10.omega( 4 ) ), 180.0, 0.001 );

		// open up psi
		movemap.set( TorsionID( 4, core::id::BB, core::id::psi_torsion ), true );
		TS_ASSERT( srfd.is_applicable( movemap, 4 ) );
		TS_ASSERT( srfd.apply( movemap, h10, 4 ) ); // psi changes here
		TS_ASSERT_EQUALS( h10.secstruct( 4 ), 'E' ); // non-intuitive behavior on part of SecstructSRFD
		TS_ASSERT_DELTA( h10.phi( 4 ), 30.0, 0.001 );
		TS_ASSERT_DELTA( h10.psi( 4 ), 40.0, 0.001 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( h10.omega( 4 ) ), 180.0, 0.001 );

		// open up omega
		movemap.set( TorsionID( 4, core::id::BB, core::id::omega_torsion ), true );
		TS_ASSERT( srfd.is_applicable( movemap, 4 ) );
		TS_ASSERT( srfd.apply( movemap, h10, 4 ) ); // omega changes here
		TS_ASSERT_EQUALS( h10.secstruct( 4 ), 'E' ); // non-intuitive behavior on part of SecstructSRFD
		TS_ASSERT_DELTA( h10.phi( 4 ), 30.0, 0.001 );
		TS_ASSERT_DELTA( h10.psi( 4 ), 40.0, 0.001 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( h10.omega( 4 ) ), 50.0, 0.001 );
	}


};
