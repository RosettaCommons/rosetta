// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/numeric_tests_that_are_easier_to_write_if_you_can_use_a_pose/cyclic_coordinate_descent.cxxtest.hh
/// @brief test suite for numeric/cyclic_coordinate_descent
/// @author Brian D. Weitzner
/// @author Jason W. Labonte

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Package headers
#include <numeric/cyclic_coordinate_descent.hh>

//Other headers
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/types.hh>

#include <numeric/wrap_angles.hh>

#include <utility/vector1.hh>

//C++ headers
#include <iostream>


class TestCyclicCoordinateDescentMath : public CxxTest::TestSuite {

private:
	core::pose::Pose pose_;

public:
	// Shared initialization goes here.
	void setUp() {
		core_init();

		using namespace core::chemical;
		using namespace core::pose;
		using core::kinematics::Edge;
		using core::kinematics::FoldTree;

		make_pose_from_sequence( pose_, "AAAAAAAAAAAAAAAAAAAAHHHH", "fa_standard" ); // spooky!

		for ( core::uint i = 1; i <= pose_.total_residue(); ++i ) {
			pose_.set_omega( i, 180.0 );
			pose_.set_phi( i, -64 );
			pose_.set_psi( i, -41 );
		}

		// Let's make a "loop" from residues 3-10 with a cutpoint at residue 6
		FoldTree ft;
		ft.add_edge( 1, 3, Edge::PEPTIDE );
		ft.add_edge( 3, 6, Edge::PEPTIDE );
		ft.add_edge( 3, 10, 1 );
		ft.add_edge( 10, 7, Edge::PEPTIDE );
		ft.add_edge( 10, pose_.total_residue(), Edge::PEPTIDE );

		// Check for stupid mistakes and set the FoldTree to the pose
		TS_ASSERT( ft.reorder( 1 ) );
		pose_.fold_tree( ft );

		// Add cutpoint variants to the loop (we use these to compute the deviation)
		core::uint cutpoint = pose_.fold_tree().cutpoint( 1 );
		add_variant_type_to_pose_residue( pose_, CUTPOINT_LOWER, cutpoint );
		add_variant_type_to_pose_residue( pose_, CUTPOINT_UPPER, cutpoint + 1 );
	}

	// Shared finalization goes here.
	void tearDown() {}

	// This helper method is reproduced from CCDLoopClosureMover
	utility::vector1< core::PointPosition >
	get_anchors( core::conformation::Residue const & residue ) const
	{
		using utility::vector1;
		using core::PointPosition;
		using namespace core::chemical;

		vector1< PointPosition > anchors;

		if ( residue.has_variant_type( CUTPOINT_LOWER ) ) {
			core::uint const last_lower_mainchain_atom( residue.mainchain_atom( residue.mainchain_atoms().size() ) );

			anchors.push_back( residue.atom( last_lower_mainchain_atom ).xyz() );
			anchors.push_back( residue.atom( "OVL1" ).xyz() );
			anchors.push_back( residue.atom( "OVL2" ).xyz() );

		} else if ( residue.has_variant_type( CUTPOINT_UPPER ) ) {
			core::uint const first_upper_mainchain_atom( residue.mainchain_atom( 1 ) );
			core::uint const second_upper_mainchain_atom( residue.mainchain_atom( 2 ) );

			anchors.push_back( residue.atom( "OVU1" ).xyz() );
			anchors.push_back( residue.atom( first_upper_mainchain_atom ).xyz() );
			anchors.push_back( residue.atom( second_upper_mainchain_atom ).xyz() );
		} else {
			TS_ASSERT( false );
		}

		return anchors;
	}

	// --------------- Test Cases --------------- //
	void test_angle_forward(){
		using numeric::ccd_angle;
		using numeric::wrap_180;
		using utility::vector1;
		using core::Angle;
		using core::DistanceSquared;
		using core::PointPosition;
		using core::pose::Pose;
		using core::Real;
		using core::Vector;

		core::uint resno = 4;
		Pose pose( pose_ );

		Angle initial_phi =  pose.phi( resno );
		core::uint cutpoint = pose.fold_tree().cutpoint( 1 );

		// Break the loop
		pose.set_phi( resno, 150 );

		// The C-side of the loop is "fixed" because resno < cutpoint
		vector1< PointPosition > moving = get_anchors( pose.residue( cutpoint ) );
		vector1< PointPosition > fixed = get_anchors( pose.residue( cutpoint + 1 ) );

		// Get the coordinates of the atom at the base of the rotation and the axis around which it rotates
		// Since we moved phi, the rotation point is about CA, and we are interested in the N->CA axis.
		PointPosition axis_atom;
		Vector axis;
		axis_atom = pose.residue( resno ).xyz( "CA" );
		axis = ( axis_atom - pose.residue( resno ).xyz( "N" ) ).normalized();

		// Awwww yeah, here we go!
		Angle alpha;
		DistanceSquared deviation;
		ccd_angle( fixed, moving, axis_atom, axis, alpha, deviation );

		// This takes care of quadrants and the like.
		Angle test_phi = wrap_180( pose.phi( resno ) + alpha );

		// If the math is correct, we should find the angle the minimizes the chain break.  Since we started with a closed
		// loop, we know that the initial value of phi is the angle that produces a chain break of 0.
		TS_ASSERT_DELTA( test_phi, initial_phi, 1E-2 );
		TS_ASSERT( deviation < 1E-6 );
	}

	void test_angle_backward(){
		using numeric::ccd_angle;
		using numeric::wrap_180;
		using utility::vector1;
		using core::Angle;
		using core::DistanceSquared;
		using core::PointPosition;
		using core::pose::Pose;
		using core::Real;
		using core::Vector;

		core::uint resno = 8;
		Pose pose( pose_ );

		Angle initial_phi =  pose.phi( resno );
		core::uint cutpoint = pose.fold_tree().cutpoint( 1 );

		// Break the loop
		pose.set_phi( resno, 150 );

		// The N-side of the loop is "fixed" because resno > cutpoint
		vector1< PointPosition > fixed = get_anchors( pose.residue( cutpoint ) );
		vector1< PointPosition > moving = get_anchors( pose.residue( cutpoint + 1 ) );

		// Get the coordinates of the atom at the base of the rotation and the axis around which it rotates
		// Since we moved phi, the rotation point is about CA, and we are interested in the N->CA axis.
		PointPosition axis_atom;
		Vector axis;
		axis_atom = pose.residue( resno ).xyz( "CA" );
		axis = ( axis_atom - pose.residue( resno ).xyz( "N" ) ).normalized();
		axis *= -1; // make sure the axis points toward the cutpoint

		// Awwww yeah, here we go!
		Angle alpha;
		DistanceSquared deviation;
		ccd_angle( fixed, moving, axis_atom, axis, alpha, deviation );

		// This takes care of quadrants and the like.
		Angle test_phi = wrap_180( pose.phi( resno ) + alpha );

		// If the math is correct, we should find the angle the minimizes the chain break.  Since we started with a closed
		// loop, we know that the initial value of phi is the angle that produces a chain break of 0.
		TS_ASSERT_DELTA( test_phi, initial_phi, 1E-2 );
		TS_ASSERT( deviation < 1E-6 );
	}

	void test_dev(){
		using numeric::ccd_angle;
		using utility::vector1;
		using core::Angle;
		using core::DistanceSquared;
		using core::PointPosition;
		using core::pose::Pose;
		using core::Real;
		using core::Vector;

		core::uint resno = 4;
		Pose pose( pose_ );

		core::uint cutpoint = pose.fold_tree().cutpoint( 1 );

		// Break the loop, but this time move several torsions
		pose.set_phi( resno, 150 );
		pose.set_psi( resno, 150 );
		pose.set_psi( resno + 4, 150 );

		// The C-side of the loop is "fixed" because resno < cutpoint
		vector1< PointPosition > moving = get_anchors( pose.residue( cutpoint ) );
		vector1< PointPosition > fixed = get_anchors( pose.residue( cutpoint + 1 ) );

		// Get the coordinates of the atom at the base of the rotation and the axis around which it rotates
		// Since we moved phi, the rotation point is about CA, and we are interested in the N->CA axis.
		PointPosition axis_atom;
		Vector axis;
		axis_atom = pose.residue( resno ).xyz( "CA" );
		axis = ( axis_atom - pose.residue( resno ).xyz( "N" ) ).normalized();

		// Awwww yeah, here we go!
		Angle alpha;
		DistanceSquared deviation;
		ccd_angle( fixed, moving, axis_atom, axis, alpha, deviation );

		// If the math is correct, we should find that the deviation of the virtual atoms in the cutpoint variants
		// matches the deviation returned by the ccd_angle function after updating the torsion
		pose.set_phi( resno, pose.phi( resno ) + alpha );

		// The ccd_angle function retuns the sum of the squared distances of the overlap atoms, so we compute that value
		// explicitly here.
		DistanceSquared sum_sq_dist( 0 );
		sum_sq_dist += pose.residue( cutpoint ).xyz( "C" ).distance_squared(
			pose.residue( cutpoint + 1 ).atom( "OVU1" ).xyz() );

		sum_sq_dist += pose.residue( cutpoint ).atom( "OVL1" ).xyz().distance_squared(
			pose.residue( cutpoint + 1 ).xyz( "N" ) );

		sum_sq_dist += pose.residue( cutpoint ).atom( "OVL2" ).xyz().distance_squared(
			pose.residue( cutpoint + 1 ).xyz( "CA" ) );

		TS_ASSERT_DELTA( deviation, sum_sq_dist, 1E-8 );
	}
};
