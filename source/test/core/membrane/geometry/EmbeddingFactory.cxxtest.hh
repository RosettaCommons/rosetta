// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/geometry/EmbeddingFactory.cxxtest.hh
///
/// @brief 		Test Suite for using the embedding factory
/// @details	Last Modified: 12/24/13
///
/// @author 	Rebecca Alford (rfalford12@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Tested Classes
#include <core/membrane/geometry/EmbeddingFactory.hh>

#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Exceptions.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/membrane/io/EmbedDefIO.hh>
#include <core/membrane/io/SpanFileIO.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>

// Resource Manager Headers
#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/LazyResourceManager.hh>
#include <basic/resource_manager/util.hh>

#include <protocols/jd2/JD2ResourceManager.hh>
#include <protocols/jd2/JD2ResourceManagerJobInputter.hh>
#include <protocols/jd2/Job.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <algorithm>
#include <cmath>

using namespace core;

/// @brief Test Suite for Membrane Embedding factory
class EmbeddingFactoryTest : public CxxTest::TestSuite {
    
public:
    
    /// @brief Setup
    void setUp()
    {
        protocols_init();
        
        // Initialize Test Data
        init_test_data();
        
        // Initialize Embedding Factories
        init_factories();
    }
    
    /// @brief teardown
    void tearDown()
    {}
    
    /// @brief Testing Init from Search and Score method
    void test_returnCorrectDefault() {
        
        using namespace core::conformation;
        
        TS_TRACE("Testing loading default embedding definition using factory ");
        factory_default_->create_and_add_embedding(true);
        
        // grab pose by rm - use util to compare data
        TS_ASSERT( compare_def_to_params( default_center_, default_normal_, default_depth_ ) );
    }
    
    /// @brief
    void test_returnCorrectPDB() {

        using namespace core::conformation;
        
        TS_TRACE("Testing loading pdb based embedding definition using factory ");
        factory_pdb_->create_and_add_embedding(true);
        
        // grab pose by rm - use util to compare data
        TS_ASSERT( compare_def_to_params( pdb_center_, pdb_normal_, pdb_depth_ ) );
        
    }
    
    /// @brief Testing correct return of topology based data
    void test_returnCorrectTopo() {
        
        using namespace core::conformation;
        
        TS_TRACE("Testing loading topology-based embedding definition using factory ");
        factory_topology_->create_and_add_embedding(true);
        
        // grab pose by rm - use util to compare data
        TS_ASSERT( compare_def_to_params( topology_center_, topology_normal_, topology_depth_ ) );
    }
    
    /// @brief testing mixed case 1
    void test_returnCorrectMix1() {
        
        using namespace core::conformation;
        
        TS_TRACE("Testing loading mixed case (type 1) embedding definition using factory ");
        factory_mix1_->create_and_add_embedding(true);
        
        // grab pose by rm - use util to compare data
        TS_ASSERT( compare_def_to_params( topology_center_, default_normal_, topology_depth_ ) );
    }
    
    /// @brief testing mixed case 2
    void test_returnCorrectMix2() {
        
        using namespace core::conformation;
        
        TS_TRACE("Testing loading mixed case (type 2) embedding definition using factory ");
        factory_mix2_->create_and_add_embedding(true);
        
        // grab pose by rm - use util to compare data
        TS_ASSERT( compare_def_to_params( pdb_center_, topology_normal_, topology_depth_ ) );
    }
    
    /// @brief testing mixed case 3
    void test_returnCorrectMix3() {
        
        using namespace core::conformation;
        
        TS_TRACE("Testing loading mixed case (type 3) embedding definition using factory ");
        factory_mix3_->create_and_add_embedding(true);
        
        // grab pose by rm - use util to compare data
        TS_ASSERT( compare_def_to_params( pdb_center_, default_normal_, topology_depth_ ) );
    }
    
    
private: // init functions
    
    /// @brief Initialize Embedding factories from resources
    void init_factories() {
        
        using namespace core::membrane::geometry;
        using namespace core::conformation::membrane;
        using namespace core::membrane::io;
        
        // Embedding Definition IO Class
        EmbedDefIO edio;
        SpanFileIO sfio;
        
        // Get pose and topology
        pose_ = new pose::Pose();
        core::import_pose::pose_from_pdb(*pose_, "core/membrane/io/1afo_test.pdb");
        SpanningTopologyOP topology = sfio.get_topology_from_spanfile("core/membrane/io/1afo_test.span");
        
        // Grabbing Embedding Data Direct from IO Classes
        EmbedConfigInfoOP def = edio.get_embedding_from_file("core/membrane/geometry/default.embed");
        EmbedConfigInfoOP pdb = edio.get_embedding_from_file("core/membrane/geometry/pdb.embed");
        EmbedConfigInfoOP by_topology = edio.get_embedding_from_file("core/membrane/geometry/topology.embed");
        
        EmbedConfigInfoOP mix1 = edio.get_embedding_from_file("core/membrane/geometry/mix1.embed");
        EmbedConfigInfoOP mix2 = edio.get_embedding_from_file("core/membrane/geometry/mix2.embed");
        EmbedConfigInfoOP mix3 = edio.get_embedding_from_file("core/membrane/geometry/mix3.embed");
        
        // Build Embedding Factories from Resources Above
        
        // Store some instances of the embedding factory :D
        factory_default_ = new EmbeddingFactory( pose_, def, topology );
        factory_pdb_ = new EmbeddingFactory( pose_, pdb, topology );
        factory_topology_ = new EmbeddingFactory( pose_, by_topology, topology );
        
        // Mixed case factories
        factory_mix1_ = new EmbeddingFactory( pose_, mix1, topology );
        factory_mix2_ = new EmbeddingFactory( pose_, mix2, topology );
        factory_mix3_ = new EmbeddingFactory( pose_, mix3, topology );
    }
    
    /// @brief Initialize Initial Conditions
    void init_test_data() {
        
        // Default Parameters
        default_center_.assign( 1, 2, 3);
        default_normal_.assign( 6, 7, 8);
        default_depth_ = 0;
        
        // PDB Based Parameters
        pdb_center_.assign( 0, 0, 0 );
        pdb_normal_.assign( 0, 0, 1 );
        pdb_depth_ = 0;
        
        // Topology Based Parameters
        // values from sc membrane test 11/24/13 - signed by RA
        topology_center_.assign( 0.98975, -0.8255, -2.32475 );
        topology_normal_.assign( -0.0904839, 0.977314, -0.191494 );
        topology_depth_ = 0;
    
    }

    /// @brief Testing helper function
    bool compare_def_to_params(
                               core::Vector center,
                               core::Vector normal,
                               core::Real depth
                               ) {
        
        using namespace core::membrane::geometry;
        
        // Check there are 81 residues in the pose
        TS_ASSERT( pose_->total_residue() == 81 );
        
        // Get residue from pose (in this test, should always be 81)
        core::conformation::Residue rsd = pose_->residue(81);
        
        // Grab Info from Embedding Definition (I think this is correct placement??)
        core::Vector const & new_center = rsd.xyz(1);
        core::Vector const & new_normal = rsd.xyz(2);
        core::Vector const & new_depth = rsd.xyz(3);
        
        // Compare Center parameters
        TS_ASSERT_DELTA( new_center.x(), center.x(), 0.0001 );
        TS_ASSERT_DELTA( new_center.y(), center.y(), 0.0001 );
        TS_ASSERT_DELTA( new_center.z(), center.z(), 0.0001 );
        
        // Compare normal parameters
        TS_ASSERT_DELTA( new_normal.x(), normal.x(), 0.0001 );
        TS_ASSERT_DELTA( new_normal.y(), normal.y(), 0.0001 );
        TS_ASSERT_DELTA( new_normal.z(), normal.z(), 0.0001 );
        
        // Compare depth parameters
        TS_ASSERT_DELTA( new_depth.y(), depth, 0.0001 );
        
        // Otherwise, criteria passed and return true
        return true;
    }
    
private: // data
    
    // Pose
    core::pose::PoseOP pose_;

    // Default Parameters
    core::Vector default_center_;
    core::Vector default_normal_;
    core::Real default_depth_;
    
    // PDB Based Parameters
    core::Vector pdb_center_;
    core::Vector pdb_normal_;
    core::Real pdb_depth_;
    
    // Topology Based Parameters
    core::Vector topology_center_;
    core::Vector topology_normal_;
    core::Real topology_depth_;
    
    // Store some instances of the embedding factory :D
    core::membrane::geometry::EmbeddingFactoryOP factory_default_;
    core::membrane::geometry::EmbeddingFactoryOP factory_pdb_;
    core::membrane::geometry::EmbeddingFactoryOP factory_topology_;
    //core::membrane::geometry::EmbeddingFactoryOP factory_score_;

    // Mixed case factories
    core::membrane::geometry::EmbeddingFactoryOP factory_mix1_;
    core::membrane::geometry::EmbeddingFactoryOP factory_mix2_;
    core::membrane::geometry::EmbeddingFactoryOP factory_mix3_;
    //core::membrane::geometry::EmbeddingFactoryOP factory_mix4_;
    
}; // class EmbeddingFactoryTest
