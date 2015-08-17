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

        // Load in pose from pdb for specific anchor point case
        positioned_pose_ = core::pose::PoseOP( new Pose() );
        pose_from_pdb( *positioned_pose_, "protocols/membrane/1C3W_TR_A.pdb" );
        
        // Load in pose from PDB containing a membrane residue with non-default position (for case 4)
        specially_positioned_pose1_ = core::pose::PoseOP( new Pose() );
        pose_from_pdb( *specially_positioned_pose1_, "protocols/membrane/2zup_model.pdb" );
        
        // Load in pose from PDB containing a membrane residue with non-default position (for case 5)
        specially_positioned_pose2_ = core::pose::PoseOP( new Pose() );
        pose_from_pdb( *specially_positioned_pose2_, "protocols/membrane/2zup_model.pdb" );
        
        // Load in a pose from PDB containing multiple membrane residues - not allowed!
        multi_mem_pose_ = core::pose::PoseOP( new Pose() );
        pose_from_pdb( *multi_mem_pose_, "protocols/membrane/1AFO_AB_multi.pdb" );
        
		// Initialize Spans from spanfile
		std::string spanfile = "protocols/membrane/1C3W_A.span";
        std::string spanfile_2zup = "protocols/membrane/2zup_model.span";

        // Testing (4) Different setups for AddMembraneMover
        //   (1) Default setup - read topology from spanfile
        //   (2) First custom setup - anchor membrane residue at an
        //       arbitrary residue in the pose and setup from a pre-defined
        //       topology object
        //   (3) Second custom setup - set a custom initial membrane position
        //   (4) Third custom setup - from existing membrane residue (points directly)
        //   (5) Fourth custom setup - find the membrane residue in Pose
		
        // (1) Setup pose with default setup
        TS_TRACE( "Setting up membrane pose with default configuration" ); 
		AddMembraneMoverOP add_memb( new AddMembraneMover( spanfile ) );
		add_memb->apply( *pose_ );
        
        // (2) Setup pose with custom anchor point & topology object
        TS_TRACE( "Setting up membrane with custom anchor point and topology object" ); 
        SpanningTopologyOP custom_topo( new SpanningTopology() );
        custom_topo->fill_from_spanfile( spanfile );
        AddMembraneMoverOP add_memb2( new AddMembraneMover( custom_topo, 50 ) );
        add_memb2->apply( *anchored_pose_ );

        // (3) Setup pose with custom membrane position
        TS_TRACE( "Setting up a pose with a custom membrane position" ); 
        Vector test_center( 10, 10, 10 ); 
        Vector test_normal( 0, 1, 0 ); // Normal along y axis
        AddMembraneMoverOP add_memb3( new AddMembraneMover( test_center, test_normal, spanfile, 0 ) ); 
		add_memb3->apply( *positioned_pose_ );
        
        // (4) Setup the pose, directly pointing to the new membrane residue
        TS_TRACE( "Setting up a pose, directly pointed to a new membrane residue already in the pose" ); 
        AddMembraneMoverOP add_memb4( new AddMembraneMover( spanfile_2zup, 334 ) );
        add_memb4->apply( *specially_positioned_pose1_ );
        
        // (5) Ask add membrane mover to go on a scavenger hunt for the membrane residue
        TS_TRACE( "Setting up a membrane pose where you need to search for the membrane residue in the pose already" ); 
        AddMembraneMoverOP add_memb5( new AddMembraneMover( spanfile_2zup, 334 ) );
        add_memb5->apply( *specially_positioned_pose2_ );
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

    /// @brief Test default membrane position setup
    void test_default_membrane_position() {

        TS_TRACE( "Test correct setup of the default membrane position: center=origin, normal along z axis" );

        // Grab current center/normal from the pose
        core::Vector current_center( pose_->conformation().membrane_info()->membrane_center() ); 
        core::Vector current_normal( pose_->conformation().membrane_info()->membrane_normal() ); 

        // Define expected center/normal
        core::Vector expected_center(0,0,0);
        core::Vector expected_normal(0,0,1);

        TS_ASSERT( position_equal_within_delta( current_center, expected_center, 0.001 ) );
        TS_ASSERT( position_equal_within_delta( current_normal, expected_normal, 0.001 ) );

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

    /// @brief Checking custom setup of initial membrane position
    void test_user_defined_membrane_position() {
	
        TS_TRACE( "Test for correct setup of a user defined membrane position" );

        // Grab current center/normal from the pose
        core::Vector current_center( positioned_pose_->conformation().membrane_info()->membrane_center() ); 
        core::Vector current_normal( positioned_pose_->conformation().membrane_info()->membrane_normal() ); 

        // Define expected center/normal
        core::Vector expected_center(10,10,10);
        core::Vector expected_normal(0,1,0);

        TS_ASSERT( position_equal_within_delta( current_center, expected_center, 0.001 ) );
        TS_ASSERT( position_equal_within_delta( current_normal, expected_normal, 0.001 ) );


    }
    
    /// @brief Checking setup from existing membrane residue (case 4 and 5)
    void test_existing_membrane_rsd() {
     
        // Existing center/normal
        Vector expected_center( -7.279,  -4.677,  -0.389 );
        Vector expected_rsd_seqpos( 334 );
        
        // Check membrane residue is at the expected position
        TS_ASSERT_EQUALS( specially_positioned_pose1_->conformation().membrane_info()->membrane_rsd_num(), expected_rsd_seqpos );
        TS_ASSERT_EQUALS( specially_positioned_pose2_->conformation().membrane_info()->membrane_rsd_num(), expected_rsd_seqpos );
        
        // Check the membrane residue has the appropriate membrane position
        position_equal_within_delta( specially_positioned_pose1_->conformation().membrane_info()->membrane_center(), expected_center, 0.001 );
        position_equal_within_delta( specially_positioned_pose2_->conformation().membrane_info()->membrane_center(), expected_center, 0.001 );
        // Not testing normal for now - test case related to a bug in relax
        
    }
    
    /// @brief Try loading in a pose with multiple membrane residues
    void test_multi_mem_pose_throws() {
        
        TS_TRACE( "Testing that an exception is thrown when trying to load in a pose with multiple membrane residues" );
        
        using namespace protocols::membrane;
        std::string spanfile = "protocols/membrane/1AFO_AB.span";
        AddMembraneMoverOP add_memb( new AddMembraneMover( spanfile ) );
        TS_ASSERT_THROWS_ANYTHING( add_memb->apply( *multi_mem_pose_ ) );
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
    core::pose::PoseOP anchored_pose_;
    core::pose::PoseOP positioned_pose_;
    core::pose::PoseOP specially_positioned_pose1_;
    core::pose::PoseOP specially_positioned_pose2_;
    core::pose::PoseOP multi_mem_pose_;
	
}; // AddMembraneMover unit test
