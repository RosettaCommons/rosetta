// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 	 core/membrane/MembraneInfo.cxxtest.hh
///
/// @brief 	 Membrane Conformation Info Object - Unit Test
/// @details The Membrane Conformation is responsible for:
///             - maintaining a correct membrane foldtree
///             - maintaining references to the membrane and embedding residues
///             - providing access to membrane related data
///
/// @note    Last Modified 3/12/14
/// @author  Rebecca Alford (rfalford12@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/conformation/membrane/MembraneInfo.hh>

// Package Headers
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/util.hh>

#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/types.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

using namespace core::kinematics;
using namespace core::conformation;
using namespace core::conformation::membrane;

class MembraneInfoTest : public CxxTest::TestSuite {
	
public: // test functions
    
    /// Test Setup Functions ////////
    
    /// @brief Setup Test
    void setUp() {
        
        using namespace core::import_pose;
        using namespace core::pose;
		
        // Initialize
        core_init();
        
        // Load in pose from pdb
        pose_ = new Pose();
        pose_from_pdb( *pose_, "core/membrane/io/1afo_test.pdb" );
        
        // Setup virtual residues
        setup_membrane_conf();
        
        // Setup the embedding residue map
        embres_map_.resize( 2 );
        embres_map_[ 1 ] = std::pair< int, int >( 1, 82 );
        embres_map_[ 2 ] = std::pair< int, int >( 41, 83 );
        
        // Setup the membrane root
        membrane_ = 81;
		
				// Construct a Membrane Info object within the conformation
				pose_->conformation().setup_membrane( embres_map_, membrane_ );

				// Extract conformation from the pose
				mp_conf_ = pose_->conformation();
        
    }
    
    /// @brief Standard Tear Down
    void tearDown() {}
	    
    ///// Test Methods /////////////
    
    /// @brief Test Membrane Fold Tree
    void test_fold_tree() {
        
        TS_TRACE("Testing membrane fold tree from membrane conformation");
        
        // Grab Fold Tree from conformation
        FoldTree const & ft = mp_conf_->fold_tree();
        
        // Check that my jump edges exist
        TS_ASSERT( ft.jump_exists( 1, 41 ) );
        TS_ASSERT( ft.jump_exists( 81, 1 ) );
        TS_ASSERT( ft.jump_exists( 1, 82 ) );
        TS_ASSERT( ft.jump_exists( 41, 83 ) );
        
        // Check that my membrane position is still the root
        TS_ASSERT( ft.is_root( 81 ) );
        
    }
    
    /// @brief Testing membrane conformation invariants
    void test_membrane_invariants() {
        
        TS_TRACE("Testing membrane conformation invariants");
        TS_ASSERT( mp_conf_->membrane()->is_membrane() );
    }
    
    /// @brief Testing membrane info access methods
    void test_membrane_info_methods() {
        
        TS_TRACE("Testing membrane info access methods");
        
        // Grabbing info from membrane positions
        core::Vector mp_center = mp_conf_->membrane()->membrane_center();
        core::Vector mp_normal = mp_conf_->membrane()->membrane_normal();
        core::Real mp_thickness = mp_conf_->membrane()->membrane_thickness();
        
        // Check center
        TS_ASSERT_EQUALS( mp_center.x(), 0 );
        TS_ASSERT_EQUALS( mp_center.y(), 0 );
        TS_ASSERT_EQUALS( mp_center.z(), 0 );
        
        // Check normal
        TS_ASSERT_EQUALS( mp_normal.x(), 0 );
        TS_ASSERT_EQUALS( mp_normal.y(), 0 );
        TS_ASSERT_EQUALS( mp_normal.z(), 1 );
        
        // Check thicnkess
        TS_ASSERT_EQUALS( mp_thickness, 30.0 );
        
    }
    
    /// @brief testing embedding info access methods
    void test_embedding_info_methods() {
        
        TS_TRACE("Testing embedding info access emthodds");
        
        // Grabbing embedding info for chain 1
        core::Vector emb_center_1 = mp_conf_->membrane()->embedding_center(1);
        core::Vector emb_normal_1 = mp_conf_->membrane()->embedding_normal(1);
        core::Real emb_depth_1 = mp_conf_->membrane()->embedding_depth(1);
        
        // Check center
        TS_ASSERT_EQUALS( emb_center_1.x(), 1 );
        TS_ASSERT_EQUALS( emb_center_1.y(), 2 );
        TS_ASSERT_EQUALS( emb_center_1.z(), 3 );
        
        // Check normal
        TS_ASSERT_EQUALS( emb_normal_1.x(), 4 );
        TS_ASSERT_EQUALS( emb_normal_1.y(), 5 );
        TS_ASSERT_EQUALS( emb_normal_1.z(), 6 );
        
        // Check thicnkess
        TS_ASSERT_EQUALS( emb_depth_1, 40.0 );
        
        // Grammbing embedding info for chain 2
        core::Vector emb_center_2 = mp_conf_->membrane()->embedding_center(2);
        core::Vector emb_normal_2 = mp_conf_->membrane()->embedding_normal(2);
        core::Real emb_depth_2 = mp_conf_->membrane()->embedding_depth(2);
        
        // Check center
        TS_ASSERT_EQUALS( emb_center_2.x(), 10 );
        TS_ASSERT_EQUALS( emb_center_2.y(), 11 );
        TS_ASSERT_EQUALS( emb_center_2.z(), 12 );
        
        // Check normal
        TS_ASSERT_EQUALS( emb_normal_2.x(), 13 );
        TS_ASSERT_EQUALS( emb_normal_2.y(), 14 );
        TS_ASSERT_EQUALS( emb_normal_2.z(), 15 );
        
        // Check thicnkess
        TS_ASSERT_EQUALS( emb_depth_2, 50.0 );
        
    }
    
    /// @brief Testing membrane conformation info getters
    void test_membrane_conf_info_accessors() {
        
        TS_TRACE("Testing membrane conformation access methods");
        
        // Check membrane info
        TS_ASSERT_EQUALS( mp_conf_->membrane()->membrane(), 81 );
        
        // Check embedding residue data
        TS_ASSERT_EQUALS( mp_conf_->membrane()->embres_map().at(1).first, 1 );
        TS_ASSERT_EQUALS( mp_conf_->membrane()->embres_map().at(2).first, 41 );
        TS_ASSERT_EQUALS( mp_conf_->membrane()->embres_map().at(1).second, 82 );
        TS_ASSERT_EQUALS( mp_conf_->membrane()->embres_map().at(2).second, 83 );
        
    }
    
    /// @brief Testing that you cannot currently insert a residue into a membrane conformation
    void test_membrane_inserton() {
        
        TS_TRACE("Testing that you cannot insert a residue into a membrane conformation");
        
        using namespace core::chemical;
        using namespace core::conformation;
        
        // Add virtual residues to pose
        // Option Setting for residue type set
        ResidueTypeSetCAP const & residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ));
        core::chemical::ResidueTypeCOPs const & rsd_type_list1( residue_set->name3_map("VRT") );
        core::chemical::ResidueType const & virtuals( *rsd_type_list1[1] );
        core::conformation::ResidueOP vrt1( core::conformation::ResidueFactory::create_residue(virtuals) );
        
       // TS_ASSERT_THROWS_ANYTHING( mp_conf_->insert_residue_by_jump( *vrt1, 40, 2, "", "", true) );
    }
    
private: // setup functions
    
    /// @brief Setup Membrane Conformaiton
    void setup_membrane_conf() {
        
        using namespace core::chemical;
        using namespace core::conformation;
        
        // Add virtual residues to pose
        // Option Setting for residue type set
        ResidueTypeSetCAP const & residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ));
        core::chemical::ResidueTypeCOPs const & rsd_type_list1( residue_set->name3_map("VRT") );
        core::chemical::ResidueType const & virtuals( *rsd_type_list1[1] );
        
        // Create three virtual residues
        core::conformation::ResidueOP vrt1( core::conformation::ResidueFactory::create_residue(virtuals) );
        core::conformation::ResidueOP vrt2( core::conformation::ResidueFactory::create_residue(virtuals) );
        core::conformation::ResidueOP vrt3( core::conformation::ResidueFactory::create_residue(virtuals) );
        
        // Setting up traditional membrane-style coordinates (vrt1 simulates the membrane)
        core::Vector center1(0, 0, 0);
        core::Vector normal1(0, 0, 1);
        core::Vector depth1( 1, 30.0, 1);
        
        vrt1->set_xyz( 2, center1 );
        vrt1->set_xyz( 1, normal1 );
        vrt1->set_xyz( 3, depth1 );
        
        // Setting up traditional embedding coordinates (vrt2, vrt3)
        core::Vector center2( 1, 2, 3);
        core::Vector normal2( 4, 5, 6);
        core::Vector depth2( 1, 40.0, 1);
        
        core::Vector center3( 10, 11, 12);
        core::Vector normal3( 13, 14, 15);
        core::Vector depth3( 1, 50.0, 1);
        
        vrt2->set_xyz( 2, center2 );
        vrt2->set_xyz( 1, normal2 );
        vrt2->set_xyz( 3, depth2 );
        
        vrt3->set_xyz( 2, center3 );
        vrt3->set_xyz( 1, normal3 );
        vrt3->set_xyz( 3, depth3 );
        
        // Add virtual residues to the existing pose
        pose_->append_residue_by_jump( *vrt1, 1, "", "", true);
        pose_->append_residue_by_jump( *vrt2, 1 );
        pose_->append_residue_by_jump( *vrt3, 41 );
        
    }
    
private: // test data
    
    // Maintain a pose
    core::pose::PoseOP pose_;
    
    // Maintain the corresponding edge list
    utility::vector1< std::pair< int, int > > embres_map_;
    
    // Maintain corresponding root
    int membrane_;
	
	// Pose Conformation
	ConformationOP mp_conf_;
};


