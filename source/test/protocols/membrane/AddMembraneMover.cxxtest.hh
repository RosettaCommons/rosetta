// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 	 protocols/membrane/AddMembraneMover.cxxtest.hh
///
/// @brief 	 Unit Test: Add membrane mover
/// @details The add membrane mover sets up an anchored membrane fold tree, is responsible for
///			 adding the membrane residue and passing correct information to MembraneInfo for setup.
///			 Last Modified: 7/8/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Package Headers
#include <protocols/membrane/AddMembraneMover.hh>

#include <core/kinematics/FoldTree.hh>

#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>

using namespace core::kinematics;
using namespace core::conformation;
using namespace core::conformation::membrane;

class AddMembraneMoverTest : public CxxTest::TestSuite {
	
public: // test functions
    
    // Test Setup Functions ///////////////////////////
    
    /// @brief Setup Test
    void setUp() {
		
		using namespace core::import_pose;
        using namespace core::pose;
		using namespace protocols::membrane;
		
        // Initialize core & options system
        core_init();
        
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
		
	}
	
	/// @brief Tear Down Test
	void tearDown() {}
	
	// Test Methods /////////////////////////////////
	
	/// @brief Check conformation invariant is true after running this mover
	void test_conformation_invariant() {
		
		TS_TRACE("Testing membrane conformation invariants");
        TS_ASSERT( pose_->conformation().is_membrane() );
		
	}
	/// @brief Check that conformation returns valid center & normal position
	void test_membrane_rsd_tracking() {
		
		TS_TRACE( "Check that membrane residue number tracks a residue of type MEM" );
		
		// Grab membrane rsd num
		core::Size resnum = pose_->conformation().membrane_info()->membrane_rsd_num();
		std::string res_type = pose_->residue( resnum ).name3();
		TS_ASSERT_EQUALS( res_type.compare( "MEM" ), 0 );
		
	}
	
	/// @brief Test residue z position method
	void test_membrane_jump_tracking() {
		
		TS_TRACE( "Check that the jump number tracks a jump containing the membrane residue" );
		
		core::Size jump = pose_->conformation().membrane_info()->membrane_jump();
		core::Size resnum = pose_->conformation().membrane_info()->membrane_rsd_num();
		core::Size expected = pose_->fold_tree().downstream_jump_residue( jump );
		TS_ASSERT_EQUALS( expected, resnum );
		
	}
	
private:
	
	core::pose::PoseOP pose_;
	
}; // AddMembraneMover unit test
