// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/geometry/geometryUtil.cxxtest.hh
///
/// @brief 		Test Suite for membrane geometry utility functions
/// @details	Last Modified: 11/14/13
///
/// @author 	Rebecca Alford (rfalford12@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Tested Classes
#include <core/membrane/geometry/util.hh>
#include <core/conformation/membrane/Exceptions.hh>

// Package Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <algorithm>
#include <cmath>

/// @brief Test Suite for Membrane Geometry Util Class
class GeometryUtilTest : public CxxTest::TestSuite {
    
public:
    
    /// @brief Setup
    void setUp()
    {
        protocols_init();
        
        // Setup residues
        create_residues();
    }
    
    /// @brief teardown
    void tearDown()
    {}
    
    /// @brief Testing Depth Calculations
    void test_mpDepth() {
        
        using namespace core::membrane::geometry;
        
        TS_TRACE("Testing membrane depth calculations per residye fron normal and center");
        
        // Using Standard Center and Normal
        core::Vector center(0, 0, 0);
        core::Vector normal(0, 0, 1);
        
        // calculating depth
        core::Real depth1 = get_mpDepth( normal, center, *rsd_zero_ );
        core::Real depth2 = get_mpDepth( normal, center, *rsd_nonzero_ );
        
        // Checking
        TS_ASSERT_DELTA( depth1, 30, 0.0001 );
        TS_ASSERT_DELTA( depth2, 34, 0.0001 );
    }

private: // methods
    
    /// @brief Create base residues
    void create_residues() {
        
        using namespace core::conformation;
        
        TS_TRACE("Creating membrane residues for testing...");
        
        bool fullatom(true);
        
        // Create ResidueTypeSet of Virtual Residues
        core::chemical::ResidueTypeSetCAP const & residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( fullatom ? core::chemical::FA_STANDARD : core::chemical::CENTROID ));
        core::chemical::ResidueTypeCOPs const & rsd_type_list( residue_set->name3_map("VRT") );
        core::chemical::ResidueType vrt_type = *rsd_type_list[1];
        
        // Create Residues for Depth tests
        rsd_zero_ = core::conformation::ResidueFactory::create_residue( vrt_type );
        rsd_nonzero_ = core::conformation::ResidueFactory::create_residue( vrt_type );
        
        core::Vector test1( 0, 0, 0);
        core::Vector test2( 2, 3, 4);
        
        rsd_zero_->set_xyz( 2, test1 );
        rsd_nonzero_->set_xyz( 2, test2 );
        
        // Setting Up Residue Checks for Equality
        core::chemical::ResidueTypeCOPs const & non_vrt_rsd_type_list( residue_set->name3_map("ALA") );
        
        // Setup Residues for .equals method test
        rsd_diff_x_ = core::conformation::ResidueFactory::create_residue( vrt_type );
        rsd_diff_y_ = core::conformation::ResidueFactory::create_residue( vrt_type );
        rsd_diff_z_ = core::conformation::ResidueFactory::create_residue( vrt_type );
        
        rsd_non_vrt_ = core::conformation::ResidueFactory::create_residue( *non_vrt_rsd_type_list[1] );
        rsd_ok_ = core::conformation::ResidueFactory::create_residue( vrt_type );
        
        // Coordinate tests
        core::Vector test3(1, 0, 0);
        core::Vector test4(0, 1, 0);
        core::Vector test5(0, 0, 1);
        
        rsd_diff_x_->set_xyz( 2, test3 );
        rsd_diff_y_->set_xyz( 2, test4 );
        rsd_diff_z_->set_xyz( 1, test5 );
        
        rsd_non_vrt_->set_xyz( 2, test1 );
        rsd_ok_->set_xyz( 2, test1 );
    }
    
private: // data
    
    // Test Residues for depth
    core::conformation::ResidueOP rsd_zero_;
    core::conformation::ResidueOP rsd_nonzero_;
    
    // Test Residues
    core::conformation::ResidueOP rsd_diff_x_;
    core::conformation::ResidueOP rsd_diff_y_;
    core::conformation::ResidueOP rsd_diff_z_;

    core::conformation::ResidueOP rsd_non_vrt_;
    
    core::conformation::ResidueOP rsd_ok_;
    
    

}; // class GeometryUtilTest
