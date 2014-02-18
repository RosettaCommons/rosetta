// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/MembraneProteinOptions.cxxtest.hh
///
/// @brief 		Test Suite for Initializing Membrane Protein Options
/// @details    CxxTest suite
///
/// @author 	Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified: 1/13/14

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Tested Classes
#include <core/membrane/MembraneProteinOptions.hh>

// Package Headers
#include <core/types.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <algorithm>

/// @brief Test Class: Membrane Protien Options (Resource Manager)
class MembraneProteinOptionsTest : public CxxTest::TestSuite {
    
public: // test methods
    
	/// @brief SetUp - Runs before each test
	void setUp()
	{
    
        using namespace core::membrane;
        
        // Initialize
		core_init();
        
        // Initilaize a new Options Object
        opts_ = new MembraneProteinOptions();
	}
    
	/// @brief tearDon - runs after each test
	void tearDown()
	{}
    
	/// Tests for Membrane Protein Options////
    
    /// @brief Test Fullatom option getter/setter pair
    void test_fullatom() {
        
        TS_TRACE("Testing membrane protein fullatom option");
        
        TS_ASSERT( opts_->fullatom() );
        opts_->fullatom(false);
        TS_ASSERT( !opts_->fullatom() );
        
    }
    
    /// @brief Test membrane protein include lips option getter/setter pair
    void test_include_lips() {
        
        TS_TRACE("Testing membrane protein include lips option");
        
        TS_ASSERT( !opts_->include_lips() );
        opts_->include_lips(true);
        TS_ASSERT( opts_->include_lips() );
    }
    
    /// @brief Test membrane protein membrane chains option getter/setter pair
    void test_membrane_chains() {
        
        TS_TRACE("Testing normal start angle");
        
        TS_ASSERT( opts_->membrane_chains().compare("") == 0 );
        opts_->membrane_chains( "chains.txt" );
        TS_ASSERT( opts_->membrane_chains().compare("chains.txt") == 0 );
        
    }
    
    private: // test data
        
        // Class to Test
        core::membrane::MembraneProteinOptionsOP opts_;
        
    }; // MembraneProteinOptionsTest

