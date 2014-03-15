// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file           core/membrane/geometry/MembraneResidueFactory.cxxtest.hh
///
/// @brief          Membrane Protein Factory Test Suite
/// @details        Test Suite for creating membrane residues/embedding residues and appending
///                 to a specified pose
///
/// @note           Last Modified: 12/1/13
/// @author         Rebecca Alford (rfalford12@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Tested Classes
#include <core/membrane/geometry/MembraneResidueFactory.hh>

// Package headers
#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/VariantType.hh>

#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <algorithm>

/// @brief Test Suite Membrane Residues Test
class MembraneResiduesTest : public CxxTest::TestSuite {
    
public:
    
    /// @brief Standard Setup Function
    void setUp() {
        
        using namespace basic::options;
        using namespace core::pose;
        using namespace core::import_pose;
        
        protocols_init();
        
        // Option - Ignore unrecognzied res (for this test pdb)
        option[ OptionKeys::in::ignore_unrecognized_res ](true);
        
        // Create an Input DataSet for residues
        center1_.assign(0, 0, 0);
        normal1_.assign(0, 0, 1);
        
        center2_.assign(1, 2, 3);
        normal2_.assign(6, 5, 4);
        
        // Load a test pose
        test_pose_ = new core::pose::Pose();
        pose_from_pdb( *test_pose_, "core/membrane/io/1afo_test.pdb");
        
    }
    
    /// @brief Testing final fold tree is valid
    void test_check_final_foldtree() {
        
        TS_TRACE("Testing new fold tree is valid");
        
        core::kinematics::FoldTree ft( test_pose_->fold_tree() );
        TS_ASSERT( ft.check_fold_tree() );
    }
    
    /// @brief Testing Center Residue Creator Method
    void test_membrane_residue()
    {
        using namespace core::pose;
        using namespace core::import_pose;
        
        TS_TRACE("Testing membrane residue definition...");
        
        // Residue found
        bool res_found = false;
        
        // Add Memrane Residues
        core::Real thickness = 30.0;
        mrf_.add_membrane_residue( center2_, normal2_, thickness, *test_pose_, true );
        
        // Going to run a loop here
        for ( core::Size i = 1; i <= test_pose_->total_residue(); ++i ) {
            
            if ( test_pose_->residue(i).type().name().compare("MEM") == 0 ) {
                
                // Check center conditions
                TS_ASSERT( test_pose_->residue(i).atom( 1 ).xyz().x() == 1 );
                TS_ASSERT( test_pose_->residue(i).atom( 1 ).xyz().y() == 2 );
                TS_ASSERT( test_pose_->residue(i).atom( 1 ).xyz().z() == 3 );
                
                // Check normal conditions
                TS_ASSERT( test_pose_->residue(i).atom( 2 ).xyz().x() == 6 );
                TS_ASSERT( test_pose_->residue(i).atom( 2 ).xyz().y() == 5 );
                TS_ASSERT( test_pose_->residue(i).atom( 2 ).xyz().z() == 4 );
                
                // Check thickness condition
                TS_ASSERT( test_pose_->residue(i).atom( 3 ).xyz().y() == 30.0 );
                
                res_found = true;
            }
        }
        
        if ( !res_found ) {
            TS_FAIL("Could not find membrane residue in the pose!");
        }
    }
    
    /// @brief testing the other one
    void test_embedding_residue()
    {
        
        using namespace core::import_pose;
        using namespace core::membrane::geometry;
        using namespace core::conformation;
        
        TS_TRACE("Testing membrane embedding definition...");
    
        bool res_found = false;
        
        // Add embedding residue by specified jump
        core::Real depth = 25.0;
        Conformation & conf = test_pose_->conformation();
        core::Size jump = conf.chain_begin(1);
        mrf_.add_embedding_residue( center1_, normal1_, depth, *test_pose_, jump, true );
        
        // Going to run the loop
        for ( core::Size i = 1; i <= test_pose_->total_residue(); ++i ) {
            
            if ( test_pose_->residue(i).type().name().compare("EMB") == 0 ) {
                
                // Check center conditions
                TS_ASSERT( test_pose_->residue(i).atom( 1 ).xyz().x() == 0 );
                TS_ASSERT( test_pose_->residue(i).atom( 1 ).xyz().y() == 0 );
                TS_ASSERT( test_pose_->residue(i).atom( 1 ).xyz().z() == 0 );
                
                // Check normal conditions
                TS_ASSERT( test_pose_->residue(i).atom( 2 ).xyz().x() == 0 );
                TS_ASSERT( test_pose_->residue(i).atom( 2 ).xyz().y() == 0 );
                TS_ASSERT( test_pose_->residue(i).atom( 2 ).xyz().z() == 1 );
                
                // Check depth condition
                TS_ASSERT( test_pose_->residue(i).atom( 3 ).xyz().y() == 25.0 );
                
                res_found = true;
            }
        }

        if ( !res_found ) {
            TS_FAIL("Could not find the embedding residue in the pose!");
        }
    }

private:
    
    // Create an instance of the MembraneResidueFactory
    core::membrane::geometry::MembraneResidueFactory mrf_;
    
    // Create a test pose
    core::pose::PoseOP test_pose_;
    
    // Storing Membrane Data (2 pairs)
    core::Vector center1_;
    core::Vector normal1_;
    
    core::Vector center2_;
    core::Vector normal2_;
    
};
