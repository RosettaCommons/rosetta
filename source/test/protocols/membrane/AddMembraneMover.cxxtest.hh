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
///			 Last Modified: 12/7/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Package Headers
#include <protocols/membrane/AddMembraneMover.hh>

#include <core/kinematics/FoldTree.hh>

#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh> 
#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/types.hh>
#include <core/conformation/membrane/types.hh>

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
        
        // Load in poses from pdb (general case)
        pose_ = core::pose::PoseOP( new Pose() );
        pose_from_pdb( *pose_, "protocols/membrane/1C3W_TR_A.pdb" );
        
        // Load in pose from pdb for specific anchor point case
        anchored_pose_ = core::pose::PoseOP( new Pose() );
        pose_from_pdb( *anchored_pose_, "protocols/membrane/1C3W_TR_A.pdb" );
        
		// Initialize Spans from spanfile
		std::string spanfile = "protocols/membrane/1C3W_A.span";
		
		// Setup membrane info object from add membrane mover
		Vector center( mem_center );
		Vector normal( mem_normal );
		
        // Setup first pose
		AddMembraneMoverOP add_memb( new AddMembraneMover( spanfile ) );
		add_memb->apply( *pose_ );
        
        // Setup second pose for anchor point - has custom topology
        SpanningTopologyOP custom_topo( new SpanningTopology() );
        custom_topo->fill_from_spanfile( spanfile );
        AddMembraneMoverOP add_memb2( new AddMembraneMover( custom_topo, 50 ) );
        add_memb2->apply( *anchored_pose_ );
		
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
    
    /// @brief Double check some initial settings from the full blown inputs
    void test_initial_spans_setup() {
        
        TS_TRACE( "Check that add membrane mover has called all of the required inputs" );
        
        // Did I load in a spanning topology?
        if ( pose_->conformation().membrane_info()->spanning_topology() != 0 ) {
            TS_ASSERT_EQUALS( pose_->conformation().membrane_info()->spanning_topology()->nspans(), 7 );
        } else {
            TS_FAIL( "Spanning topology not initialized. Abort!" );
        }
        
        // Check that I did not include a lipsfile quite yet
        TS_TRACE( "Check that I have not yet initialized a lipophilicity object yet" );
        TS_ASSERT( !pose_->conformation().membrane_info()->include_lips() );
    }
	
	/// @brief Test that the membrane jump number is at the expected position
    /// during initialization
	void test_membrane_jump_tracking() {
		
		TS_TRACE( "Check that the jump number tracks a jump containing the membrane residue" );
		
		core::Size jump = pose_->conformation().membrane_info()->membrane_jump();
        core::Size expected( 1 );
		TS_ASSERT_EQUALS( expected, jump );
		
	}
    
    /// @brief Check the initial configuration of the foldtre contains
    /// a single membrane residue connected by jump to the protein and that
    /// it is the root
    void test_initial_foldtree() {
        
        TS_TRACE( "Check that the membrane foldtree is setup contianing a membrane residue at the root of a simple tree" );
        
        // Redundant double checking - ensure the membrane residue is the last
        // residue in the initial pose
        core::Size mprsd = pose_->conformation().membrane_info()->membrane_rsd_num();
        core::Size expected_resnum( 223 );
        TS_TRACE("Check that the memrbane residue is located at the end of the pose");
        TS_ASSERT_EQUALS( expected_resnum, mprsd );
        TS_ASSERT_EQUALS( pose_->total_residue(), expected_resnum );
        
        // Check that the root of the pose is the membrane residue num
        core::Size expected_root( 223 );
        core::Size given_root( pose_->fold_tree().root() );
        TS_TRACE("Check that the root of the foldtree is the membrane residue");
        TS_ASSERT_EQUALS( given_root, expected_root );
        
        // Check that the membrane residue is connected to the first
        // residue in the protein
        core::Size jump = pose_->conformation().membrane_info()->membrane_jump();
        core::Size expected_upstream( 223 );
        core::Size given_upstream(  pose_->fold_tree().upstream_jump_residue( jump ) );
        core::Size expected_downstream( 1 );
        core::Size given_downstream( pose_->fold_tree().downstream_jump_residue( jump ) );
        
        // Check that the upstream and downstream resnums match
        TS_TRACE( "Checking upstream (root) residue numbers match in the membrane jump");
        TS_ASSERT_EQUALS( given_upstream, expected_upstream );
        TS_TRACE( "Checking downstream (pose first residue) residue number matches in the mmebrane jump" );
        TS_ASSERT_EQUALS( given_downstream, expected_downstream );
        
    }
    
    /// @brief Test for custom anchored pose setup (anchor the membrane residue to
    /// some abitrary point on the pose
    void test_anchored_foldtree() {
        
        TS_TRACE( "Check that the anchored foldtree (custom) has a correct setup to start" );
        
        // Check that the root of the pose is the membrane residue num
        core::Size expected_root( 223 );
        core::Size given_root( anchored_pose_->fold_tree().root() );
        TS_TRACE("Check that the root of the foldtree is the membrane residue for anchored foldtree");
        TS_ASSERT_EQUALS( given_root, expected_root );
        
        // Check that the membrane residue is connected to residue 50
        // in the protein
        core::Size jump = anchored_pose_->conformation().membrane_info()->membrane_jump();
        core::Size expected_upstream( 223 );
        core::Size given_upstream(  anchored_pose_->fold_tree().upstream_jump_residue( jump ) );
        core::Size expected_downstream( 50 );
        core::Size given_downstream( anchored_pose_->fold_tree().downstream_jump_residue( jump ) );
        
        // Check that the upstream and downstream resnums match
        TS_TRACE( "Checking upstream (root) residue numbers match in the membrane jump");
        TS_ASSERT_EQUALS( given_upstream, expected_upstream );
        TS_TRACE( "Checking downstream (pose first residue) residue number matches in the mmebrane jump" );
        TS_ASSERT_EQUALS( given_downstream, expected_downstream );
    }
	
private:
	
	core::pose::PoseOP pose_;
    core::pose::PoseOP anchored_pose_;
	
}; // AddMembraneMover unit test
