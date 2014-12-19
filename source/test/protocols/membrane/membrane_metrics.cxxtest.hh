// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 	 protocols/membrane/metrics.cxxtest.hh
///
/// @brief 	 Unit Test: Membrane model analysis metrics
/// @details Methods in this file are used for assessment of membrane protein models. Currently
///          includes rmsd mehtods over the transmembrane spanning regions.
///			 Last Modified: 10/24/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Package Headers
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/metrics.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/types.hh>
#include <core/conformation/membrane/types.hh>

// Utility Headers
#include <utility/vector1.hh>

/// @brief Unit tests for membrane RMSD methods
class MembraneMetricsTest : public CxxTest::TestSuite {
    
public: // test functions
    
    // Test Setup Functions ///////////////////////////
    
    /// @brief Setup Test
    void setUp() {
        
        using namespace core::import_pose;
        using namespace core::pose;
		using namespace core::conformation::membrane;
        using namespace protocols::membrane;
		
        // Initialize core & options system
        core_init();
        
        // Load in native pose from pdb
        native_pose_ = core::pose::PoseOP( new Pose() );
        pose_from_pdb( *native_pose_, "protocols/membrane/1afo_in.pdb" );
        
        // Load in test pose from pdb
        test_pose_ = core::pose::PoseOP( new Pose() );
        pose_from_pdb( *test_pose_, "protocols/membrane/1afo_decoy.pdb" );
        
        // Initialize Spans from spanfile
        std::string spanfile = "protocols/membrane/1afo_tr.span";
        
        // Setup membrane info object from add membrane mover
        Vector center( mem_center );
        Vector normal( mem_normal );
        
        AddMembraneMoverOP add_memb( new AddMembraneMover( center, normal, spanfile, 1 ) );
        add_memb->apply( *native_pose_ );
        add_memb->apply( *test_pose_ );
        
    }
    
    /// @brief Tear Down Test
    void tearDown() {}
    
    // Test Methods /////////////////////////////////
    
    /// @brief Calculate membrane backbone rmsd with superposition
    void test_membrane_bb_rmsd_with_super() {
        
        using namespace protocols::membrane;
        
        TS_TRACE( "Calculating membrane backbone rmsd with superposition" );
        
        core::Real rms = mem_bb_rmsd_with_super( *native_pose_, *test_pose_ );
        TS_ASSERT_DELTA( rms, 0.0005, 0.0003 );
    
    }
    
    /// @brief Calculate membrane backbone rmsd with superposition
    void test_membrane_bb_rmsd_no_super() {
        
        using namespace protocols::membrane;
        
        TS_TRACE( "Calculating membrane backbone rmsd without superposition" );
        
        core::Real rms = mem_bb_rmsd_no_super( *native_pose_, *test_pose_ );
        TS_ASSERT_DELTA( rms, 3.8042, 0.0003 );
        
    }

    /// @brief Calculate membrane all atom rmsd with superposition
    void test_membrane_bb_rmsd_with_super_allatom() {
        
        using namespace protocols::membrane;
        
        TS_TRACE( "Calculating membrane allatom rmsd with superposition" );
        
        core::Real rms = mem_all_atom_rmsd_with_super( *native_pose_, *test_pose_ );
        TS_ASSERT_DELTA( rms, 0.5052, 0.0003 );
        
    }
    
    /// @brief Calculate membrane all atom rmsd without superposition
    void test_membrane_bb_rmsd_no_super_allatom() {
        
        using namespace protocols::membrane;
        
        TS_TRACE( "Calculating membrane backbone rmsd without superposition" );
        
        core::Real rms = mem_all_atom_rmsd_no_super( *native_pose_, *test_pose_ );
        TS_ASSERT_DELTA( rms, 3.9376, 0.0003 );
        
    }

    
private:
    
    core::pose::PoseOP native_pose_;
    core::pose::PoseOP test_pose_;
    
}; // Membrane metrics unit test
