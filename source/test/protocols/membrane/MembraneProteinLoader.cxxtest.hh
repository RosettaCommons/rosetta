// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		protocols/membrane/MembraneProteinLoader.cxxtest.hh
///
/// @brief 		Test Suite for group of classes that load in an embedding definition resource
/// @details	Tests span file loader functionality using generic loader class
///
/// @author     Rebecca Alford


// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/membrane/MembraneMover.hh>

// Project Headers
#include <core/membrane/MembraneProteinLoader.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

/// @brief Test Suite for Membrane Embedding factory
class MembraneProteinLoaderTest : public CxxTest::TestSuite {
    
public:
    
    
    /// @brief Setup
    void setUp()
    {
        using namespace basic::options;
        protocols_init();
        
        // Register Appropriate Options
        option[OptionKeys::in::membrane](true);
        option[OptionKeys::in::file::membrane_chains]("protocols/membrane/chains.txt");
        option[OptionKeys::in::file::fullatom](true);
        
        // set resource definition files here
        option[OptionKeys::jd2::resource_definition_files]("protocols/membrane/membrane_full.xml");
        option[OptionKeys::out::overwrite](true);
        
        try {
            
            using namespace protocols::membrane;
            using namespace protocols::jd2;
            using namespace core::import_pose;
            
            // Initialize Membrane Mover
            MembraneMoverOP mp = new MembraneMover();
            mp->set_full(false);
            JobDistributor::get_instance()->go(mp);
            
            // Get Membrane Pose from the Membrane Mover
            membrane_protein_ = mp->get_pose();
            
        } catch ( utility::excn::EXCN_Base const & e ) {
            std::cout << "caught exception " << e.msg() << std::endl;
        }
        
    }
    
    /// @brief teardown
    void tearDown()
    {}
    
    /// @brief Initial Test Code
    void test_membrane_protein() {
        
        TS_TRACE("Placeholder");
        TS_ASSERT( true );
    }
    
private: // data
    
    // Resulting Membrane Protein
    core::pose::PoseOP membrane_protein_;
    
}; // class EmbeddingFactoryTest

