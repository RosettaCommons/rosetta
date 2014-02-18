// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/scoring/MembranePenalties.cxxtest.hh
///
/// @brief 		Test Suite for calculating membrane penalties
/// @details	Tests using topology based data
///
/// @author 	Rebecca Alford (rfalford12@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Tested Classes
#include <core/membrane/scoring/MembranePenalties.hh>

#include <core/membrane/util/definitions.hh>
#include <core/membrane/util/Exceptions.hh>

// Package Headers
#include <core/types.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <algorithm>

/// @brief Test Class Embedding Definition Data Loader
class MembranePenaltiesTests :  public CxxTest::TestSuite {
public:

    /// @brief Setup Function - Runs before every test
    void setUp()
    {
        // Initialize
        core_init();
        
        // Initializing Embedding COnfig Info Objects
        init_ = new core::membrane::util::EmbedConfigInfo();
        final_ = new core::membrane::util::EmbedConfigInfo();
        
        // Construct Membrane Penalties object
        //penalty_ = new core::membrane::scoring::MembranePenalties();
        
        // Initialize Starting data
        init_->normal.x() = -0.0904839;
        init_->normal.y() = 0.977314;
        init_->normal.z() = -0.191494;
        
        init_->center.x() = 0.98975;
        init_->center.y() = -0.8255;
        init_->center.z() = -2.32475;
        
        // Initialize Final data
        final_->normal.x() = -0.141;
        final_->normal.y() = 0.872;
        final_->normal.z() = -0.468;
        
        final_->center.x() = -1.005;
        final_->center.y() = 3.879;
        final_->center.z() = -1.558;
    }
    
    /// @brief Tear Down Function - runs after every test
    void tearDown()
    {}
        
    /// @brief Testing Starting TMH Penalty
    void test_start_tmh_penalty() {
        
        // Create score and ref score
        core::Real score = 0;
        core::Real & tm_proj ( score );
        
        // Pass by reference
        core::Vector const & normal( init_->normal );
        core::Vector const & center( init_->center );
        
        // Calculate penalty
        penalty_.tm_projection_penalty("core/membrane/io/1afo", "test", normal, center, tm_proj );
        
        // Check
        TS_ASSERT_DELTA( 200, tm_proj, 0.0001);
    }
    
    /// @brief Testing final TMH Penalty
    void test_final_tmh_penalty() {
        
        // Create score and ref score
        core::Real score = 0;
        core::Real & tm_proj ( score );
        
        // Pass by reference
        core::Vector const & normal( final_->normal );
        core::Vector const & center( final_->center );
        
        // Calculate penalty
        penalty_.tm_projection_penalty("core/membrane/io/1afo", "test", normal, center, tm_proj );
        
        // Check
        TS_ASSERT_DELTA( 100, tm_proj, 0.0001);
        
    }
    /// @brief Testing starting non helix in membrane penalty
    void test_start_nonhelix_penalty() {
        
        // Create score and ref score
        core::Real score = 0;
        core::Real & non_helix ( score );
        
        // Pass by reference
        core::Vector const & normal( init_->normal );
        core::Vector const & center( init_->center );
        
        // Calculate penalty
        penalty_.non_helix_in_membrane_penalty("core/membrane/io/1afo", "test", normal, center, non_helix );
        
        // Check
        TS_ASSERT_DELTA( 420, non_helix, 0.0001);
    }
    
    /// @brief Testing final non helix in membrane penalty
    void test_final_nonhelix_penalty() {
        
        // Create score and ref score
        core::Real score = 0;
        core::Real & non_helix ( score );
        
        // Pass by reference
        core::Vector const & normal( final_->normal );
        core::Vector const & center( final_->center );
        
        // Calculate penalty
        penalty_.non_helix_in_membrane_penalty("core/membrane/io/1afo", "test", normal, center, non_helix );
        
        // Check
        TS_ASSERT_DELTA( 370, non_helix, 0.0001);
        
    }
    
    /// @brief Testing starting termini penalty
    void test_start_termini_penalty() {
        
        // Create score and ref score
        core::Real score = 0;
        core::Real & termini ( score );
        
        // Pass by reference
        core::Vector const & normal( init_->normal );
        core::Vector const & center( init_->center );
        
        // Calculate penalty
        penalty_.termini_penalty("core/membrane/io/1afo", "test", normal, center, termini );
        
        // Check
        TS_ASSERT_DELTA( 150, termini, 0.0001);
    }
    
    /// @brief Testing final terimi penalty
    void test_final_termini_penalty() {
        
        // Create score and ref score
        core::Real score = 0;
        core::Real & termini ( score );
        
        // Pass by reference
        core::Vector const & normal( final_->normal );
        core::Vector const & center( final_->center );
        
        // Calculate penalty
        penalty_.termini_penalty("core/membrane/io/1afo", "test", normal, center, termini );
        
        // Check
        TS_ASSERT_DELTA( 150, termini, 0.0001);
        
    }
    
private:
    
    // Starting and Final Embedding Data
    core::membrane::util::EmbedConfigInfoOP init_;
    core::membrane::util::EmbedConfigInfoOP final_;
    
    // Membrane Penalties Object
    core::membrane::scoring::MembranePenalties penalty_;

};
