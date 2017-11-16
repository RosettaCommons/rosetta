// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/membrane/symmetry/SymmetricAddMembraneMover.cxxtest.hh
///
/// @brief      Unit test for symmetric add membrane mover
/// @details Checks that the membrane residue is anchored by the appropriate
///             virtual residue, is the root of the foldtree, number of virtuals
///             is appropriately updated, jump in the foldtree matches the
///             expected setup and that pose is still symmetric & membrane after
///             runtime. Since we are only supporting C symmetries and the
///             general add membrane restrictions are already tested in the
///             AddMembraneMover unit test, only checking this on 1afo for now.
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Package Headers
#include <protocols/membrane/symmetry/SymmetricAddMembraneMover.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>

#include <protocols/membrane/AddMembraneMover.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/conformation/membrane/MembraneInfo.hh>

#include <core/pose/symmetry/util.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.membrane.symmetry.SymmetricAddMembraneMover.cxxtest");

using namespace core::kinematics;
using namespace core::conformation;
using namespace core::conformation::membrane;
using namespace core::conformation::symmetry;
using namespace protocols::simple_moves::symmetry;

class SymmetricAddMembraneMoverTest : public CxxTest::TestSuite {

public: // test functions

	// Test Setup Functions ///////////////////////////

	/// @brief Setup Test
	void setUp() {

		using namespace core::import_pose;
		using namespace core::pose;
		using namespace protocols::membrane::symmetry;
		using namespace protocols::simple_moves;
		using namespace core::conformation::symmetry;

		// Initialize core & options system
		core_init();

		// Load in poses from pdb (general case)
		symmetric_pose_ = core::pose::PoseOP( new Pose() );
		pose_from_file( *symmetric_pose_, "protocols/membrane/symmetry/1afo_tr_input.pdb" , core::import_pose::PDB_file);

		// Load symmdata object from symmetry definition file
		SymmDataOP symm_data = SymmDataOP( new SymmData() );
		symm_data->read_symmetry_data_from_file( "protocols/membrane/symmetry/1afo_tr.c2.symm" );

		// Create asymmetric pose from symmetry data object
		SetupForSymmetryMoverOP setup_for_symm = SetupForSymmetryMoverOP( new SetupForSymmetryMover( symm_data ) );
		setup_for_symm->apply( *symmetric_pose_ );

		// Store path to spanfile and do symmetric add membrane
		std::string spanfile = "protocols/membrane/symmetry/1afo_tr.span";
		SymmetricAddMembraneMoverOP symm_add = SymmetricAddMembraneMoverOP( new SymmetricAddMembraneMover( spanfile ) );
		symm_add->apply( *symmetric_pose_ );

	}

	/// @brief Tear Down Test
	void tearDown() {}

	// Test Methods /////////////////////////////////

	/// @brief Check postconditions: pose is still symmetric AND a membrane protein
	void test_symm_add_postconditions() {

		TR <<  "Testing that after SymmetricAddMembrane operation, pose is still symmetric AND a membrane protein"  << std::endl;
		TS_ASSERT( symmetric_pose_->conformation().is_membrane() );
		TS_ASSERT( core::pose::symmetry::is_symmetric( *symmetric_pose_ ) );

	}

	/// @brief Check that the conformation contains the correct # of virtual residues
	/// and preserves the number of subunits and number of residues in the monomer
	/// are preserved
	void test_preserved_symmetry() {

		using namespace core::conformation::symmetry;

		TR <<  "Check that the conformation still contains the correct # of virtual residues, preserves the number of subunits, and number of residues in the monomeric unit"  << std::endl;

		SymmetricConformation & symm_conf ( dynamic_cast< SymmetricConformation & > ( symmetric_pose_->conformation()) );

		// 1afo setup with C2 symmetry
		TS_ASSERT_EQUALS( symm_conf.Symmetry_Info()->subunits(), 2 );

		// Expected nres subunut = 40 (2 chains of 1afo, 40 res each)
		TS_ASSERT_EQUALS( symm_conf.Symmetry_Info()->get_nres_subunit(), 40 );

		// Expected Num virtuals: Base + 2x V0 + 2* V1 + MEM = 6
		TS_ASSERT_EQUALS( symm_conf.Symmetry_Info()->num_virtuals(), 6 );
	}

	/// @brief Check that the membrane virtual is appropriately anchored to the
	/// V0,1 base residue and is the root of the foldtree
	void test_appropriate_foldtree() {

		TR <<  "Check that the membrane is still the root of the foldtree and appropriately attached to the framework of symmetry virtual residues"  << std::endl;

		// Redundant double checking - ensure the membrane residue is the last
		// residue in the initial pose
		core::Size mprsd = symmetric_pose_->conformation().membrane_info()->membrane_rsd_num();
		core::Size expected_resnum( 86 );
		TR << "Check that the memrbane residue is located at the end of the pose (total residue includes symmetry virtuals here!)" << std::endl;
		TS_ASSERT_EQUALS( expected_resnum, mprsd );
		TS_ASSERT_EQUALS( symmetric_pose_->size(), expected_resnum );

		// Check that the root of the pose is the membrane residue num
		core::Size expected_root( 86 );
		core::Size given_root( symmetric_pose_->fold_tree().root() );
		TR << "Check that the root of the foldtree is the membrane residue" << std::endl;
		TS_ASSERT_EQUALS( given_root, expected_root );

		// Check that the membrane residue is connected to the first
		// residue in the protein
		core::Size jump = symmetric_pose_->conformation().membrane_info()->membrane_jump();
		core::Size expected_upstream( 86 );
		core::Size given_upstream( symmetric_pose_->fold_tree().upstream_jump_residue( jump ) );
		core::Size expected_downstream( 81 );
		core::Size given_downstream( symmetric_pose_->fold_tree().downstream_jump_residue( jump ) );

		// Check that the upstream and downstream resnums match
		TR <<  "Checking upstream (root) residue numbers match in the membrane jump" << std::endl;
		TS_ASSERT_EQUALS( given_upstream, expected_upstream );
		TR <<  "Checking downstream (pose first residue) residue number matches in the mmebrane jump"  << std::endl;
		TS_ASSERT_EQUALS( given_downstream, expected_downstream );

	}


	/// @brief Test that the membrane jump number is at the expected position
	/// during initialization
	void test_membrane_jump_tracking() {

		TR <<  "Check that the jump number tracks a jump containing the membrane residue"  << std::endl;

		core::Size jump = symmetric_pose_->conformation().membrane_info()->membrane_jump();
		core::Size expected( 7 );
		TS_ASSERT_EQUALS( expected, jump );

	}

private:

	core::pose::PoseOP symmetric_pose_;

}; // SymmetricAddMembraneMover unit test
