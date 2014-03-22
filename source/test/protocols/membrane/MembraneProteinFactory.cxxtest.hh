// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 	 protocols/memrane/MembraneProteinFactory.cxxtest.hh
///
/// @brief 	 MembraneProteinFactory Test Suite
/// @details The membrane protein factory creates a single pose from various membrane proteins
///			 loaded on the front end and initialized as membrane proteins. This single framework
///			 will then be passed off to the MembraneHub (which coordinates I/O) and sent back to the protocol it was
///			 called from (usually in pose loading)
///
/// @note    Last Modified 1/4/14
/// @author  Rebecca Alford (rfalford12@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/membrane/MembraneUnitTestMover.hh>

// Project Headers
#include <core/conformation/membrane/SpanningTopology.hh>

#include <core/conformation/membrane/definitions.hh>
#include <core/conformation/membrane/Exceptions.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
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

using namespace core;

/// @brief Test Suite for Membrane Embedding factory
class MembraneProteinFactoryTest : public CxxTest::TestSuite {
    
public:
    
    // force rebuild
    
    /// @brief Setup
    void setUp()
    {
        using namespace basic::options;
        
        protocols_init();
        
        // Dummy string for something else
        option[ OptionKeys::in::file::membrane_chains ]("protocols/membrane/chains.txt");
        
        // set resource definition files here
        option[OptionKeys::jd2::resource_definition_files]("protocols/membrane/membrane.xml");
        option[OptionKeys::out::overwrite](true);
        
        try {
            
            using namespace protocols::membrane;
            using namespace protocols::jd2;
            using namespace core::import_pose;
            
            // Initialize Membrane Mover
            MembraneUnitTestMoverOP mp = new MembraneUnitTestMover();
            mp->register_options();
            mp->init_from_cmd();
            JobDistributor::get_instance()->go(mp);
            
            // Get Membrane Pose from the Membrane Mover
            membrane_protein_ = mp->get_membrane_pose();
            
        } catch ( utility::excn::EXCN_Base const & e ) {
            std::cout << "caught exception " << e.msg() << std::endl;
        }

    }
    
    /// @brief teardown
    void tearDown()
    {}
    
    /// @brief Testing pose appending
    void test_membrane_protein() {
        
        using namespace core::conformation::membrane;
        using namespace core::kinematics;
        using namespace core::membrane;
        
        TS_TRACE("Testing memrbane pose construction from multiple chains");
        
        // Check number of chains
        core::conformation::Conformation conf = membrane_protein_->conformation();
        TS_ASSERT_EQUALS( conf.num_chains(), 3 );
        
        TS_TRACE("Testing Membrane Virtual residue construction/placement");
        
        // Check Embedding + Membrane Residu
        TS_ASSERT( membrane_protein_->residue(81).type().name().compare("MEM") == 0 );
        TS_ASSERT( membrane_protein_->residue(82).type().name().compare("EMB") == 0 );
        TS_ASSERT( membrane_protein_->residue(83).type().name().compare("EMB") == 0 );
        
        TS_TRACE("Testing membrane protein constains an accessible membrane conformation object");
        
        // Cast conformation to membrane conformaiton
		TS_ASSERT( membrane_protein_->conformation().membrane()->is_membrane() ); // check dst invariants
        
        TS_TRACE("Test base case of membrane fold tree construction");
        FoldTree ft( membrane_protein_->fold_tree() );
        TS_ASSERT( ft.check_fold_tree() );
        TS_ASSERT( ft.is_root( 81 ) );
        
        TS_TRACE("Test number of chains in membrane spanning topology");
        utility::vector1< SpanningTopology > topology = membrane_protein_->conformation().membrane()->spanning_topology();
        TS_ASSERT_EQUALS( topology.size(), 2 );
        TS_ASSERT_EQUALS( topology[ 1 ].total_residue_in_span_file(), 40 );
        TS_ASSERT_EQUALS( topology[ 2 ].total_residue_in_span_file(), 40 );
    }
    
private: // data
    
    // Resulting Membrane Protein
    core::pose::PoseOP membrane_protein_;
    
}; // class MembraneProteinFactoryTest

