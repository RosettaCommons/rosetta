// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 	 core/conformation/membrane/MembraneConformation.cxxtest.hh
///
/// @brief 	 Unit Test: Methods in Conformation for Membrane Proteins
/// @details Check that methods in conformation relating to membrane proteins are working properly
///			 This includes the following methods:
///				(1) #membrane_center()
///				(2) #membrane_normal()
///				(3) membrane fold tree setup (is reasonable?)
///				(4) Computing relative z position of residues and atoms
///			 Last Modified: 7/8/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Package Headers
#include <protocols/membrane/AddMembraneMover.hh>

#include <core/kinematics/FoldTree.hh>

#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/membrane_new.OptionKeys.gen.hh>

using namespace core::kinematics;
using namespace core::conformation;
using namespace core::conformation::membrane;

class MembraneConformationTest : public CxxTest::TestSuite {
	
public: // test functions
    
    // Test Setup Functions ///////////////////////////
    
    /// @brief Setup Test
    void setUp() {
		
		using namespace basic::options;
		using namespace core::import_pose;
        using namespace core::pose;
		using namespace protocols::membrane;
		
        // Initialize core & options system
        core_init();
		
		// Make sure view in pymol is turned off
		option[ OptionKeys::membrane_new::view_in_pymol ](false);
        
        // Load in pose from pdb
        pose_ = new Pose();
        pose_from_pdb( *pose_, "protocols/membrane/1C3W_TR_A.pdb" );
        
		// Initialize Spans from spanfile
		std::string spanfile = "protocols/membrane/1C3W_A.span";
	
		// Setup membrane info object from add membrane mover
		Vector center( 0, 0, 0 );
		Vector normal( 0, 0, 1 );
		
		AddMembraneMoverOP add_memb = new AddMembraneMover( center, normal, spanfile, 1 );
		add_memb->apply( *pose_ );
		
		// Show the current setup
		pose_->conformation().show_membrane();
	}
	
	/// @brief Tear Down Test
	void tearDown() {}
	
	// Test Methods /////////////////////////////////
	
	/// @brief Check that conformation returns valid center & normal position
	void test_normal_center() {
		
		TS_TRACE( "Check that conformation returns a reasonable normal and center" );
		
		Vector center = pose_->conformation().membrane_center();
		Vector expected_center( 0, 0, 0 );
		TS_ASSERT( position_equal_within_delta( center, expected_center, 0.0001 ) );
		
		Vector normal = pose_->conformation().membrane_normal();
		Vector expected_normal( 0, 0, 1 );
		TS_ASSERT( position_equal_within_delta( normal, expected_normal, 0.0001 ) );
		
	}
	
	/// @brief Test residue z position method
	void test_residue_z_position() {
		
		TS_TRACE( "Check computing relative residue z position in the membrane" );
		
		core::Real expected_z( -22.089 );
		core::Real z( pose_->conformation().residue_z_position( 1 ) );
		TS_ASSERT_DELTA( z, expected_z, 0.0001 );
	}
	
	/// @brief Test atom z position method
	void test_atom_z_position() {
	
		TS_TRACE( "Check computing relative atom z position in the membrane" );
		
		core::Real expected_z( -22.089 );
		core::Real z( pose_->conformation().atom_z_position( 1, 2 ) );
		TS_ASSERT_DELTA( z, expected_z, 0.0001 );
	}
	
	/// @brief Check the fold tree has a reasonable setup
	void test_reasonable_foldtree() {
		
		TS_TRACE( "Checking that the defualt setup for the membrane fold tree is a reasonable fold tree" );
		
		// Grab fold tree from the pose & check
		FoldTree ft = pose_->fold_tree();
		ft.show( std::cout );
		TS_ASSERT( pose_->conformation().membrane()->check_membrane_fold_tree( ft ) );
	}
	
	/// @brief Check the fold tree is unreasonble if membrane rsd is not a jump point
	void test_unreasonable_foldtree() {
		
		TS_TRACE( "Checking that foldtree invariant returns false on a fold tree where membrane residue is not the jump point" );
		
		// Make new simple tree
		FoldTreeOP ft = new FoldTree();
		ft->simple_tree( pose_->total_residue() );
		TS_ASSERT(! pose_->conformation().membrane()->check_membrane_fold_tree( *ft ) );
		
	}
	
	/// @brief Position equal within delta (helper method)
	bool position_equal_within_delta( Vector a, Vector b, Real delta ) {
		
		TS_ASSERT_DELTA( a.x(), b.x(), delta );
		TS_ASSERT_DELTA( a.y(), b.y(), delta );
		TS_ASSERT_DELTA( a.z(), b.z(), delta );
		
		return true;
	}
	
	
private:
	
	core::pose::PoseOP pose_;
	
}; // MembraneConformation unit test
