// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 	 protocols/relax/membrane/MPFastRelaxFoldTree.cxxtest.hh
///
/// @brief 	 Unit Test: Foldtree setup for the memrbane fast relax mover
/// @details Checks that the mpfast relax mover sets up an appropriate initial foldtree:
///          membrane residue attached to the protein/molecule COM, COM is the root, membrane
///          moveable, all existing jumps between chains preserved. Tests a single chain case
///          and multiple chain case whose COM occurs in a given chain.
///			 Last Modified: 12/7/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Package Headers
#include <protocols/relax/membrane/MPFastRelaxMover.hh>
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
using namespace protocols::relax::membrane;

class MPFastRelaxFoldTreeTest : public CxxTest::TestSuite {
    
public: // test functions
    
    // Test Setup Functions ///////////////////////////
    
    /// @brief Setup Test
    void setUp() {
        
        using namespace core::import_pose;
        using namespace core::pose;
        using namespace protocols::membrane;
        
        // Initialize core & options system
        core_init();
        
        // Load in 1C3W (Bacteriorhodopsin) Test case: single chain
        pose_1c3w_ = core::pose::PoseOP( new Pose() );
        pose_from_pdb( *pose_1c3w_, "protocols/membrane/1C3W_TR_A.pdb" );
        std::string spanfile = "protocols/membrane/1C3W_A.span";
        Vector center( 0, 0, 0 );
        Vector normal( 0, 0, 1 );
        AddMembraneMoverOP add_memb1( new AddMembraneMover( center, normal, spanfile ) );
        add_memb1->apply( *pose_1c3w_ );
        
        // Load in 2mpn (Inner membrane protein YgaP) Test case: two chains (docking...?)
        pose_2mpn_ = core::pose::PoseOP( new Pose() );
        pose_from_pdb( *pose_2mpn_, "protocols/relax/membrane/2mpn_tr_native.pdb" );
        std::string spanfile2 = "protocols/relax/membrane/2mpn_tr.span";
        AddMembraneMoverOP add_memb2( new AddMembraneMover( center, normal, spanfile2 ) );
        add_memb2->apply( *pose_2mpn_ );
        
        // TODO: Add symmetric, asymmetric test cases. Will also get expanded if there
        // is a version of symmetric membrane relax in the works.
        
        // Initialize the fast relax mover
        relax_mover_ = MPFastRelaxMoverOP( new MPFastRelaxMover() );
    }
    
    /// @brief Tear Down Test
    void tearDown() {}
    
    // Test Methods /////////////////////////////////
    
    /// @brief Check appropriate setup of a single chain foldtree
    /// @details Tests the 1C3W case
    void test_single_chain_com_anchored_foldtree() {
        
        relax_mover_->setup_relax_foldtree( *pose_1c3w_ );
     
        TS_TRACE( "Check the resulting foldtree for a single chain pose" );
        pose_1c3w_->fold_tree().show( std::cout );
        
        TS_TRACE( "Check that the root is the pose COM and the memrbane residue is the downstream jump" );
        
        // Check root
        core::Size expected_root( 86 );
        core::Size given_root( pose_1c3w_->fold_tree().root() );
        TS_ASSERT_EQUALS( expected_root, given_root );
        
        // Check jump
        core::Size membrane_jump( pose_1c3w_->conformation().membrane_info()->membrane_jump() );
        core::Size expected_jump( 1 );
        TS_ASSERT_EQUALS( expected_jump, membrane_jump );
        
        // Check downstream position
        core::Size downstream_pos( pose_1c3w_->fold_tree().downstream_jump_residue( membrane_jump ) );
        core::Size expected_downstream( 223 );
        TS_ASSERT_EQUALS( expected_downstream, downstream_pos );
        
        // Check num jumps in the foldtree is appropriately
        core::Size num_jumps( pose_1c3w_->fold_tree().num_jump() );
        core::Size expected_num_jumps( 1 );
        TS_ASSERT_EQUALS( expected_num_jumps, num_jumps );
        
    }
    
    /// @brief Check appropriate setup of a single chain foldtree
    /// @details Tests the 1C3W case
    void test_double_chain_com_anchored_foldtree() {
        
        relax_mover_->setup_relax_foldtree( *pose_2mpn_ );
        
        TS_TRACE( "Check the resulting foldtree for a double chain pose" );
        pose_2mpn_->fold_tree().show( std::cout );
      
        // Check root
        core::Size expected_root( 84 );
        core::Size given_root( pose_2mpn_->fold_tree().root() );
        TS_ASSERT_EQUALS( expected_root, given_root );
        
        // Check jump
        core::Size membrane_jump( pose_2mpn_->conformation().membrane_info()->membrane_jump() );
        core::Size expected_jump( 1 );
        TS_ASSERT_EQUALS( expected_jump, membrane_jump );
        
        // Check downstream position
        core::Size downstream_pos( pose_2mpn_->fold_tree().downstream_jump_residue( membrane_jump ) );
        core::Size expected_downstream( 137 );
        TS_ASSERT_EQUALS( expected_downstream, downstream_pos );
        
        // Check num jumps in the foldtree is appropriately
        core::Size num_jumps( pose_2mpn_->fold_tree().num_jump() );
        core::Size expected_num_jumps( 2 );
        TS_ASSERT_EQUALS( expected_num_jumps, num_jumps );
    }
    
    
private:
    
    // Store mp fast relax mover
    MPFastRelaxMoverOP relax_mover_;
    
    // Store poses
    core::pose::PoseOP pose_2mpn_;
    core::pose::PoseOP pose_1c3w_;
    
}; // MPFastRelaxFoldTree unit test
