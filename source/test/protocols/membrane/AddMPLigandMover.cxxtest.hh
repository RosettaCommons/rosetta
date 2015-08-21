// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/membrane/AddMPLigandMoverTest.cxxtest.hh
///
/// @brief  Unit test for a ligand to a membrane framework refinement case
/// @details  Accommodate membrane protein ligand in the membrane framework by
///    reorganizing the current foldtree. Resulting foldtree will
///    keep the membrane attached to the COM and ligand to the closest
///    binding pocket residue, provided in the constructor.
///
/// @author  Rebecca Faye Alford (rfalford12@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/membrane/AddMPLigandMover.hh>

// Package Headers
#include <protocols/membrane/AddMembraneMover.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>

#include <core/kinematics/FoldTree.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>

using namespace core;
using namespace core::kinematics;
using namespace core::conformation;
using namespace core::conformation::membrane;

/// @brief Unit test for add MP ligand mover
class AddMPLigandMoverTest : public CxxTest::TestSuite {

public: // test functions

	// Test Setup Functions ///////////////////////////

	/// @brief Setup Test
	void setUp() {

		using namespace core::import_pose;
		using namespace core::pose;
		using namespace protocols::membrane;

		// Initialize core & options system
		core_init_with_additional_options("-in:file:extra_res_fa protocols/membrane/RET.params");

		// Load in poses from pdb (general case)
		pose_with_ligand_ = core::pose::PoseOP( new Pose() );
		pose_from_pdb( *pose_with_ligand_, "protocols/membrane/3PXO_A_tr.pdb" );

		// Initialize Spans from spanfile
		std::string spanfile = "protocols/membrane/3PXO_A_tr.span";

		// Add Membrane Residue to Pose
		AddMembraneMoverOP add_memb( new AddMembraneMover( spanfile ) );
		add_memb->apply( *pose_with_ligand_ );

		// Add Ligand mover
		/// 118 = closest residue (by inspection in PyMOL), pose.total_rsd() = RET position
		AddMPLigandMoverOP add_ligand( new AddMPLigandMover( 118, pose_with_ligand_->total_residue()-1 ) );
		add_ligand->apply( *pose_with_ligand_ );

	}

	/// @brief Tear Down Test
	void tearDown() {}

	// Test Methods /////////////////////////////////

	/// @brief Check the protein is still a membrane protein after this mvoer runs
	void test_membrane_invariant() {

		TS_TRACE("Testing membrane conformation invariants");
		TS_ASSERT( pose_with_ligand_->conformation().is_membrane() );

	}

	/// @brief Check the protein is still a membrane protein after this mvoer runs
	void test_membrane_plus_ligand_foldtree() {

		TS_TRACE( "Check correct setup of the membrane foldtree including the ligand" );

		// Check that the root of the pose is still the COM
		core::Size expected_root( 121 );
		core::Size given_root( pose_with_ligand_->fold_tree().root() );
		TS_TRACE("Check that the root of the foldtree is the protein center of mass");
		TS_ASSERT_EQUALS( given_root, expected_root );

		// Check that the membrane residue is connected to the protein center of mass
		core::Size jump = pose_with_ligand_->conformation().membrane_info()->membrane_jump();
		core::Size expected_mp_upstm( 121 );
		core::Size given_mp_upstm(  pose_with_ligand_->fold_tree().upstream_jump_residue( jump ) );
		core::Size expected_mp_dwnstm( 328 );
		core::Size given_mp_dwnstm( pose_with_ligand_->fold_tree().downstream_jump_residue( jump ) );

		// Check that the upstream and downstream resnums match for the membrane jump
		TS_TRACE( "Checking upstream (root) residue numbers match in the membrane jump");
		TS_ASSERT_EQUALS( given_mp_upstm, expected_mp_upstm );
		TS_TRACE( "Checking downstream (pose first residue) residue number matches in the mmebrane jump" );
		TS_ASSERT_EQUALS( given_mp_dwnstm, expected_mp_dwnstm );\

			// Check that the membrane residue is connected to the protein center of mass
			core::Size ligand_jump( 2 );
		core::Size expected_ligand_upstm( 327 );
		core::Size given_ligand_upstm(  pose_with_ligand_->fold_tree().upstream_jump_residue( ligand_jump ) );
		core::Size expected_ligand_dwnstm( 118 );
		core::Size given_ligand_dwnstm( pose_with_ligand_->fold_tree().downstream_jump_residue( ligand_jump ) );

		// Check that the upstream and downstream resnums match for the ligand jump
		TS_TRACE( "Checking upstream (ligand closest rsd) residue numbers match in the ligand jump");
		TS_ASSERT_EQUALS( expected_ligand_upstm, given_ligand_upstm );
		TS_TRACE( "Checking downstream (ligand) residue number matches in the ligand jump" );
		TS_ASSERT_EQUALS( expected_ligand_dwnstm, given_ligand_dwnstm );
	}

	/// @brief Position equal within delta (helper method)
	bool position_equal_within_delta( Vector a, Vector b, Real delta ) {

		TS_ASSERT_DELTA( a.x(), b.x(), delta );
		TS_ASSERT_DELTA( a.y(), b.y(), delta );
		TS_ASSERT_DELTA( a.z(), b.z(), delta );

		return true;
	}

private:

	core::pose::PoseOP pose_with_ligand_;

}; // AddMembraneMover unit test



