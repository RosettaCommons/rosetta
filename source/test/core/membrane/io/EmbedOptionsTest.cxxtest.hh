// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/io/EmbedOptionsTest.cxxtest.hh
///
/// @brief 		Test Suite for Initializing Embedding Search Parameters Object
/// @details    CxxTest suite
///
/// @author 	Rebecca Alford (rfalford12@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Tested Classes
#include <core/membrane/io/EmbedSearchParamsOptions.hh>
#include <core/conformation/membrane/definitions.hh>
#include <core/conformation/membrane/definitions_util.hh>

// Package Headers
#include <core/types.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <algorithm>

/// @brief Test Class: Span File loader
class EmbedSearchParamsOptionsTests : public CxxTest::TestSuite {

public: // test methods
    
	/// @brief SetUp - Runs before each test
	void setUp()
	{
        // Initialize
		core_init();
	}
    
	/// @brief tearDon - runs after each test
	void tearDown()
	{}
    
	/// Tests for IO Class ////
    
    /// @brief Test Normal Search Getter/setter pair
    void test_normal_search() {
        
        TS_TRACE("Testing normal search");
        
        TS_ASSERT( !params.normal_search() );
        params_.set_normal_search(true);
        TS_ASSERT( params_.normal_search() );
        
    }
    
    /// @brief Test normal delta angle getter/setter pair with bounds
    void test_normal_delta_angle() {
        
        TS_TRACE("Testing normal delta angle");
        
        TS_ASSERT_EQUALS( params_.normal_delta_angle(), 0);
        
        // Try to set out of bounds
        TS_ASSERT_THROWS_ANYTHING( params_.set_normal_delta_angle(-1) );
        TS_ASSERT_THROWS_ANYTHING( params_.set_normal_delta_angle(370) );
        
        // Try to set out of bounds
        params_.set_normal_delta_angle(10);
        TS_ASSERT_EQUALS( params_.normal_delta_angle(), 10 );
    }
    
    /// @brief Test normal start angle with getter/setter pair and bounds
    void test_normal_delta_angle() {
        
        TS_TRACE("Testing normal start angle");
        
        TS_ASSERT_EQUALS( params_.normal_start_angle(), 0);
        
        // Try to set out of bounds
        TS_ASSERT_THROWS_ANYTHING( params_.set_normal_start_angle(-1) );
        TS_ASSERT_THROWS_ANYTHING( params_.set_normal_start_angle(370) );
        
        // Try to set out of bounds
        params_.set_normal_start_angle(10);
        TS_ASSERT_EQUALS( params_.normal_start_angle(), 10 );
        
    }
    
    /// @brief Test normal max angle with getter/setter pair and bounds
    void test_normal_max_angle() {
        
        TS_TRACE("Testing normal max angle");
        
        TS_ASSERT_EQUALS( params_.normal_max_angle(), 0);
        
        // Try to set out of bounds
        TS_ASSERT_THROWS_ANYTHING( params_.set_normal_max_angle(-1) );
        TS_ASSERT_THROWS_ANYTHING( params_.set_normal_max_angle(370) );
        
        // Try to set out of bounds
        params_.set_normal_max_angle(10);
        TS_ASSERT_EQUALS( params_.normal_max_angle(), 10 );
    }
    
    /// @brief center search getter/setter pair
    void test_center_search() {
        
        TS_TRACE("Testing center search");
        TS_ASSERT( !params_.center_search() );
        params_.set_center_search(true);
        TS_ASSERT( params_.center_search() );
    
    /// @brief Test center max delta with bounds
    void test_center_max_delta() {
        
        TS_TRACE("Testing center max delta");
        TS_ASSERT_EQUALS( params_.center_max_delta(), 0 );
        
        // Try to set out of bounds
        TS_ASSERT_THROWS_ANYTHING( params_.set_normal_max_angle(-1) );
        
        // Try to set out of bounds
        params_.set_center_max_delta(10);
        TS_ASSERT_EQUALS( params_.center_max_delta(), 10 );
    }
    
    /// @brief Test center magnitude
    void test_center_mag() {
        
        TS_TRACE("Testing center mag");
        TS_ASSERT_EQUALS( params_.center_mag(), 0 );
        
        // Try to set out of bounds
        params_.set_center_mag(10);
        TS_ASSERT_EQUALS( params_.center_mag(), 10 );
        
    }
    
    /// @brief Test normal magnitude
    void test_normal_mag() {
        
        TS_ASSERT("Testing normal mag");
        TS_ASSERT_EQUALS( params_.normal_mag(), 0 );
        
        // Try to set out of bounds
        params_.set_normal_mag(10);
        TS_ASSERT_EQUALS( params_.normal_mag(), 10 );
        
    }
    
    /// @brief Test normal cycles
    void test_normal_cycles() {
        
        TS_TRACE("Testing normal cycles");
        TS_ASSERT_EQUALS( params_.normal_cycles(), 0 );
        
        // Try to set out of bounds
        TS_ASSERT_THROWS_ANYTHING( params_.set_normal_cycles(-1) );
        
        // Try to set out of bounds
        params_.set_normal_cycles(10);
        TS_ASSERT_EQUALS( params_.normal_cycles(), 10 );
        
    }
    
    /// @brief Testing penalties
    void test_penalties() {
        
        TS_TRACE("Testing penalties");
        TS_ASSERT( !params_.penalties() );
        params_.set_penalties(true);
        TS_ASSERT( params_.penalties() );
        
    }
    
    /// @brief Testing No interpolate mpair
    void test_no_interpolate_mpair() {
        
        TS_TRACE("Testing no interpolate mpair");
        TS_ASSERT( !params_.no_interpolate_mpair() );
        params_.set_no_interpolate_mpair(true);
        TS_ASSERT( params_.no_interpolate_mpair() );
        
    }
    
private: // test data
    
	// Class to Test
    core::membrane::io::EmbedSearchParamsOptions params_;
    
};