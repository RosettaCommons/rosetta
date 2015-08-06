// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 	 protocols/membrane/MembraneUtil.cxxtest.hh
/// @brief 	 Unit Test: Protocols-level RosettaMP Framework calculations
///			 Last Modified: 6/10/15
/// @author  Rebecca Alford (rfalford12@gmail.com)
/// @author  JKLeman (julia.koehler1982@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Package Headers
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/util.hh>

#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/geometry/Embedding.hh>

#include <core/conformation/Conformation.hh> 
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh>

#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/pose/Pose.hh> 
#include <core/import_pose/import_pose.hh>
#include <core/types.hh> 

#include <core/conformation/membrane/Exceptions.hh>

// Utility Headers
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

using namespace core;
using namespace core::conformation;
using namespace core::conformation::membrane;
using namespace core::pose;
using namespace core::kinematics;
using namespace protocols::membrane;
using namespace protocols::membrane::geometry;


/// @brief Unit Test suite for protocols-level membrane util class
class MembraneUtil : public CxxTest::TestSuite {
    
public: // test functions
    
    // Test Setup Functions ///////////////////////////
    
    /// @brief Setup the unit test
    void setUp() {
        
        using namespace core::import_pose;
    
        // Initialize core & options system
        core_init();
        
        // Test Cases for calc angles & moveable membranes:
        // 1. TM domain of the M2 proton channel (single helix)
        m2_pose_ = PoseOP( new Pose() );
        pose_from_pdb( *m2_pose_, "protocols/membrane/1mp6_transformed.pdb" );
        AddMembraneMoverOP add_memb1 = AddMembraneMoverOP( new AddMembraneMover( "protocols/membrane/1mp6.span" ) );
        add_memb1->apply( *m2_pose_ );

        // 2. Glycophorin A (two helices, oriented 'somewhat' opposite one another)
        glpA_pose_ = PoseOP( new Pose() );
        pose_from_pdb( *glpA_pose_, "protocols/membrane/1AFO_AB.pdb" );
        AddMembraneMoverOP add_memb2 = AddMembraneMoverOP( new AddMembraneMover( "protocols/membrane/1AFO_AB.span" ) );
        add_memb2->apply( *glpA_pose_ );
		
		// Test Cases for membrane rmsd calculations
		// 3. Native Glycophprin A
		native_pose_ = core::pose::PoseOP( new Pose() );
        pose_from_pdb( *native_pose_, "protocols/membrane/1afo_in.pdb" );
		
		// 4. Transformed glycophorin A
        test_pose_ = core::pose::PoseOP( new Pose() );
        pose_from_pdb( *test_pose_, "protocols/membrane/1afo_decoy.pdb" );
		
        AddMembraneMoverOP add_memb3( new AddMembraneMover( "protocols/membrane/1afo_tr.span", 1 ) );
        add_memb3->apply( *native_pose_ );
        add_memb3->apply( *test_pose_ );
		
    }

    /// @brief Tear down the unit test
    void tearDown() {}
    
    // Test calc helix axis method
    void test_calc_helix_axis() {
        
        TS_TRACE( "Testing the helix axes are correctly calculated for each span in the test set" );
		
		// Calculate helix axis for each transembrane span
		Vector m2_first_span( calc_helix_axis( *m2_pose_, 1 ) );
		Vector glpA_first_span( calc_helix_axis( *glpA_pose_, 1 ) );
		Vector glpA_second_span( calc_helix_axis( *glpA_pose_, 2 ) );
		
		Vector m2_first_expected( -3.208, 1.1903, 30.47033 );
        Vector glpA_first_expected( 1.409, -3.2603, 23.7866 );
		Vector glpA_second_expected( -0.886, 13.6927, 23.1553 );
		
		position_equal_within_delta( m2_first_span, m2_first_expected, 0.005 );
		position_equal_within_delta( glpA_first_span, glpA_first_expected, 0.005 );
		position_equal_within_delta( glpA_second_span, glpA_second_expected, 0.005 );

    }
    
    
    // Test Calc Angles Method
    
    /// @brief Calculate angles from single helix pose and two helix pose
    void test_calc_helix_angles() {
        
        TS_TRACE( "Testing method for calculating the angle between a single helix and the membrane normal in a single helix pose" );
		
		// Calculate helix tilt angles & compare
        Real m2_first_angle( calc_helix_tilt_angle( *m2_pose_, 1 ) );
		Real glpA_first_angle( calc_helix_tilt_angle( *glpA_pose_, 1 ) );
		Real glpA_second_angle( calc_helix_tilt_angle( *glpA_pose_, 2 ) );
		
		TS_ASSERT_DELTA( m2_first_angle, 6.41, 0.005 );
		TS_ASSERT_DELTA( glpA_first_angle, 8.49, 0.005 );
		TS_ASSERT_DELTA( glpA_second_angle, 30.65, 0.005 );
		
    }
	
    
    // Test Angle RMSD method
    
    /// @brief Calculate RMSD between reference and new angle
    void test_angle_rmsd_method() {
        
        TS_TRACE( "Testing method for calcuating the rms between a measured and reference angle value" );
        
        // Test ref = 0, measured = 0
        Real rms1( calc_angle_rmsd(0, 0) );
        TS_ASSERT_EQUALS( rms1, 0.0 );
        
        // Test measured = 0, ref = positive
        Real rms2( calc_angle_rmsd( 0, 5 ) );
        TS_ASSERT_DELTA( rms2, 3.535, 0.005 );
        
        // Test measured = 0, ref = negative
        Real rms3( calc_angle_rmsd( 0, -5 ) );
        TS_ASSERT_DELTA( rms3, 3.535, 0.005 );
        
        // Test normal - both positive
        Real rms4( calc_angle_rmsd( 3, 5 ) );
        TS_ASSERT_DELTA( rms4, 1.414, 0.005 );
        
        // Test normal - both negative
        Real rms5( calc_angle_rmsd( -3, -5 ) );
        TS_ASSERT_DELTA( rms5, 1.414, 0.005 );
        
    }
    
    ///// These tests are separate and should be moved down into membrane info after debugging ////
    // There are more test cases here, but I'm only including the ones in use for now. Will extend later
    void test_is_fixed_on_moveable_memb() {
        
        TS_TRACE( "Test is_membrane_fixed method on a moveable membrane" );
        
		// Set up new foldtree for m2 rooted at the previous downstream residue
		// of the membrane jump
        Size m2_jump( m2_pose_->conformation().membrane_info()->membrane_jump() );
        Size m2_downstream( m2_pose_->conformation().fold_tree().downstream_jump_residue( m2_jump ) );
        FoldTreeOP m2_foldtree = FoldTreeOP( new FoldTree( m2_pose_->conformation().fold_tree() ) );
        m2_foldtree->reorder( m2_downstream );
		m2_pose_->fold_tree( *m2_foldtree );
        m2_foldtree->show( std::cout );
		TS_ASSERT( !is_membrane_fixed( *m2_pose_ ) );
		
		// Set up new foldtree for glpA rooted at the previous downstream residue
		// of the membrane jump
		Size glpA_jump( glpA_pose_->conformation().membrane_info()->membrane_jump() );
        Size glpA_downstream( glpA_pose_->conformation().fold_tree().downstream_jump_residue( glpA_jump ) );
		FoldTreeOP glpA_foldtree = FoldTreeOP( new FoldTree( glpA_pose_->conformation().fold_tree() ) );
		glpA_foldtree->reorder( glpA_downstream );
		glpA_pose_->fold_tree( *glpA_foldtree );
		glpA_foldtree->show( std::cout );
		TS_ASSERT( !is_membrane_fixed( *glpA_pose_ ) );
        
    }
    
    void test_is_fixed_on_fixed_memb() {
        
        TS_TRACE( "Test is_fixed method on a fixed membrane" );

        // When the pose is loaded in via AddMembraneMover, by default the foldtree
		// enforces a fixed membrane position
        TS_ASSERT( is_membrane_fixed( *m2_pose_ ) );
        TS_ASSERT( is_membrane_fixed( *glpA_pose_ ) );
        
    }
    
    void test_is_independently_moveable_on_fixed_memb() {
        
        TS_TRACE( "Test is_independently_moveable on a fixed membrane" );
		
		TS_ASSERT( !is_membrane_moveable_by_itself( *m2_pose_ ) );
		TS_ASSERT( !is_membrane_moveable_by_itself( *glpA_pose_ ) );

    }
    
    void test_is_independently_moveable_on_independently_moveable_memb() {
        
        TS_TRACE( "Test is_independently_moveable on a moveable membrane" );

		// Set up new foldtree for m2 rooted at the previous downstream residue
		// of the membrane jump
        Size m2_jump( m2_pose_->conformation().membrane_info()->membrane_jump() );
        Size m2_downstream( m2_pose_->conformation().fold_tree().downstream_jump_residue( m2_jump ) );
        FoldTreeOP m2_foldtree = FoldTreeOP( new FoldTree( m2_pose_->conformation().fold_tree() ) );
        m2_foldtree->reorder( m2_downstream );
		m2_pose_->fold_tree( *m2_foldtree );
        m2_foldtree->show( std::cout );
		TS_ASSERT( is_membrane_moveable_by_itself( *m2_pose_ ) );
		
		// Set up new foldtree for glpA rooted at the previous downstream residue
		// of the membrane jump
		Size glpA_jump( glpA_pose_->conformation().membrane_info()->membrane_jump() );
        Size glpA_downstream( glpA_pose_->conformation().fold_tree().downstream_jump_residue( glpA_jump ) );
		FoldTreeOP glpA_foldtree = FoldTreeOP( new FoldTree( glpA_pose_->conformation().fold_tree() ) );
		glpA_foldtree->reorder( glpA_downstream );
		glpA_pose_->fold_tree( *glpA_foldtree );
		glpA_foldtree->show( std::cout );
		TS_ASSERT( is_membrane_moveable_by_itself( *glpA_pose_ ) ); 
        
    }
	
	/// @brief Calculate membrane backbone rmsd with superposition
    void test_membrane_bb_rmsd_with_super() {
		
        TS_TRACE( "Calculating membrane backbone rmsd with superposition" );
        core::Real rms = mem_bb_rmsd_with_super( *native_pose_, *test_pose_ );
        TS_ASSERT_DELTA( rms, 0.0005, 0.0003 );
    
    }
    
    /// @brief Calculate membrane backbone rmsd with superposition
    void test_membrane_bb_rmsd_no_super() {
		
        TS_TRACE( "Calculating membrane backbone rmsd without superposition" );
        core::Real rms = mem_bb_rmsd_no_super( *native_pose_, *test_pose_ );
        TS_ASSERT_DELTA( rms, 3.8042, 0.0003 );
        
    }

    /// @brief Calculate membrane all atom rmsd with superposition
    void test_membrane_bb_rmsd_with_super_allatom() {
		
        TS_TRACE( "Calculating membrane allatom rmsd with superposition" );
        core::Real rms = mem_all_atom_rmsd_with_super( *native_pose_, *test_pose_ );
        TS_ASSERT_DELTA( rms, 0.5052, 0.0003 );
        
    }
    
    /// @brief Calculate membrane all atom rmsd without superposition
    void test_membrane_bb_rmsd_no_super_allatom() {
		
        TS_TRACE( "Calculating membrane backbone rmsd without superposition" );
        core::Real rms = mem_all_atom_rmsd_no_super( *native_pose_, *test_pose_ );
        TS_ASSERT_DELTA( rms, 3.9376, 0.0003 );
        
    }

    // compute_structure_based_embedding
    void test_compute_structure_based_embedding1() {
        
        TS_TRACE("Test compute_structure_based_embedding 1");
        using namespace protocols::membrane;
        using namespace protocols::membrane::geometry;
        
        // 1AFO
        // read in pose and create topology object
        TS_TRACE("1AFO");
        Pose pose1;
        core::import_pose::pose_from_pdb( pose1, "protocols/membrane/geometry/1AFO_.pdb" );

        // create membrane pose
        AddMembraneMoverOP addmem1( new AddMembraneMover( "protocols/membrane/geometry/1AFO__tr.span" ) );
        addmem1->apply( pose1 );

        // define vectors and object
        Vector center1(0.53225, 0.361, 0.095);
        Vector normal1(-0.849116, -0.512025, 0.129741);
        
        // compute embedding
        EmbeddingDefOP embed1( compute_structure_based_embedding( pose1 ) );
        
        // compare
        TS_ASSERT( position_equal_within_delta( embed1->center(), center1, 0.001 ) );
        TS_ASSERT( position_equal_within_delta( embed1->normal(), normal1, 0.001 ) );
        

    }

    // compute_structure_based_embedding
    void test_compute_structure_based_embedding2() {

        TS_TRACE("Test compute_structure_based_embedding 2");
        using namespace protocols::membrane;
        using namespace protocols::membrane::geometry;
        
        // 1BL8
        TS_TRACE("1BL8");
        Vector center( 0, 0, 0 );
        Vector normal( 0, 0, 1 );
        Pose pose2;
        core::import_pose::pose_from_pdb( pose2, "protocols/membrane/geometry/1BL8_.pdb" );
        AddMembraneMoverOP addmem2( new AddMembraneMover( "protocols/membrane/geometry/1BL8__tr.span" ) );
        addmem2->apply(pose2);
        Vector center2(73.9421, 26.7549, 24.4493);
        Vector normal2(0.384026, -0.0403822, 0.922439);
        EmbeddingDefOP embed2( compute_structure_based_embedding( pose2 ) );
        TS_ASSERT( position_equal_within_delta( embed2->center(), center2, 0.001 ) );
        TS_ASSERT( position_equal_within_delta( embed2->normal(), normal2, 0.001 ) );

    }

    // compute_structure_based_embedding
    void test_compute_structure_based_embedding3() {
        
        TS_TRACE("Test compute_structure_based_embedding 3");
        using namespace protocols::membrane;
        using namespace protocols::membrane::geometry;

        // 1QJP - beta-barrel
        TS_TRACE("1QJP");
        Vector center( 0, 0, 0 );
        Vector normal( 0, 0, 1 );
        Pose pose3;
        core::import_pose::pose_from_pdb( pose3, "protocols/membrane/geometry/1QJP_.pdb" );
        AddMembraneMoverOP addmem3( new AddMembraneMover( "protocols/membrane/geometry/1QJP__tr.span" ) );
        addmem3->apply(pose3);
        Vector center3(31.2161, 16.9685, 37.6119);
        Vector normal3(0.877926, -0.47167, 0.0822946);
        EmbeddingDefOP embed3( compute_structure_based_embedding( pose3 ) );
        TS_ASSERT( position_equal_within_delta( embed3->center(), center3, 0.001 ) );
        TS_ASSERT( position_equal_within_delta( embed3->normal(), normal3, 0.001 ) );
        
    }

    // compute_structure_based_embedding
    void test_compute_structure_based_embedding4() {
        
        TS_TRACE("Test compute_structure_based_embedding 4");
        using namespace protocols::membrane;
        using namespace protocols::membrane::geometry;

        // 2BS2
        TS_TRACE("2BS2");
        Vector center( 0, 0, 0 );
        Vector normal( 0, 0, 1 );
        Pose pose4;
        core::import_pose::pose_from_pdb( pose4, "protocols/membrane/geometry/2BS2_CF.pdb" );
        AddMembraneMoverOP addmem4( new AddMembraneMover( "protocols/membrane/geometry/2BS2_CF_tr.span" ) );
        addmem4->apply(pose4);
        Vector center4(21.4326, 6.0464, -41.0573);
        Vector normal4(0.0060, 0.0117348, 0.9999);
        EmbeddingDefOP embed4( compute_structure_based_embedding( pose4 ) );
        TS_ASSERT( position_equal_within_delta( embed4->center(), center4, 0.001 ) );
        TS_ASSERT( position_equal_within_delta( embed4->normal(), normal4, 0.001 ) );
        
    }

    // compute_structure_based_embedding
    void test_compute_structure_based_embedding5() {
        
        TS_TRACE("Test compute_structure_based_embedding 5");
        using namespace protocols::membrane;
        using namespace protocols::membrane::geometry;
        
        // 2MPN
        TS_TRACE("2MPN");
        Vector center( 0, 0, 0 );
        Vector normal( 0, 0, 1 );
        Pose pose5;
        core::import_pose::pose_from_pdb( pose5, "protocols/membrane/geometry/2MPN_.pdb" );
        AddMembraneMoverOP addmem5( new AddMembraneMover( "protocols/membrane/geometry/2MPN__tr.span" ) );
        addmem5->apply(pose5);
        Vector center5(0.3645, 3.66025, 41.345);
        Vector normal5(0.006956, 0.996657, 0.0813998);
        EmbeddingDefOP embed5( compute_structure_based_embedding( pose5 ) );
        TS_ASSERT( position_equal_within_delta( embed5->center(), center5, 0.001 ) );
        TS_ASSERT( position_equal_within_delta( embed5->normal(), normal5, 0.001 ) );
        
    }

    // compute_structure_based_embedding
    void test_compute_structure_based_embedding6() {
        
        TS_TRACE("Test compute_structure_based_embedding 6");
        using namespace protocols::membrane;
        using namespace protocols::membrane::geometry;
        
        // 2OAR
        TS_TRACE("2OAR");
        Vector center( 0, 0, 0 );
        Vector normal( 0, 0, 1 );
        Pose pose6;
        core::import_pose::pose_from_pdb( pose6, "protocols/membrane/geometry/2OAR_.pdb" );
        AddMembraneMoverOP addmem6( new AddMembraneMover( "protocols/membrane/geometry/2OAR__tr.span" ) );
        addmem6->apply(pose6);
        Vector center6(18.8453, 122.117, 1.079);
        Vector normal6(-0.755094, -0.655578, 0.00710981);
        EmbeddingDefOP embed6( compute_structure_based_embedding( pose6 ) );
        TS_ASSERT( position_equal_within_delta( embed6->center(), center6, 0.001 ) );
        TS_ASSERT( position_equal_within_delta( embed6->normal(), normal6, 0.001 ) );
        
    }

    // compute_structure_based_embedding
    void test_compute_structure_based_embedding7() {
        
        TS_TRACE("Test compute_structure_based_embedding 7");
        using namespace protocols::membrane;
        using namespace protocols::membrane::geometry;
        
        // 2UUH
        TS_TRACE("2UUH");
        Vector center( 0, 0, 0 );
        Vector normal( 0, 0, 1 );
        Pose pose7;
        core::import_pose::pose_from_pdb( pose7, "protocols/membrane/geometry/2UUH__tr.pdb" );
        AddMembraneMoverOP addmem7( new AddMembraneMover( "protocols/membrane/geometry/2UUH__tr.span" ) );
        addmem7->apply(pose7);
        Vector center7(-0.000166667, -0.000125, 0.295625);
        Vector normal7(0, 0, 1);
        EmbeddingDefOP embed7( compute_structure_based_embedding( pose7 ) );
        TS_ASSERT( position_equal_within_delta( embed7->center(), center7, 0.001 ) );
        TS_ASSERT( position_equal_within_delta( embed7->normal(), normal7, 0.001 ) );
        
        
    }

    // compute_structure_based_embedding
    void test_compute_structure_based_embedding8() {
        
        TS_TRACE("Test compute_structure_based_embedding 8");
        using namespace protocols::membrane;
        using namespace protocols::membrane::geometry;

        // 3PXO
        TS_TRACE("3PXO");
        Vector center( 0, 0, 0 );
        Vector normal( 0, 0, 1 );
        Pose pose8;
        core::import_pose::pose_from_pdb( pose8, "protocols/membrane/geometry/3PXO_.pdb" );
        AddMembraneMoverOP addmem8( new AddMembraneMover( "protocols/membrane/geometry/3PXO__tr.span" ) );
        addmem8->apply(pose8);
        Vector center8(-36.1201, -7.59636, 37.6713);
        Vector normal8(-0.9862, -0.164797, 0.0158378);
        EmbeddingDefOP embed8( compute_structure_based_embedding( pose8 ) );
        TS_ASSERT( position_equal_within_delta( embed8->center(), center8, 0.001 ) );
        TS_ASSERT( position_equal_within_delta( embed8->normal(), normal8, 0.001 ) );
        
    }

    // check vector for reasonable size
    void test_check_vector() {
        
        TS_TRACE("Test check vector");
                
        // define vectors and object
        Vector v1(1, 2, 1000);
        Vector v2(4, 5000, 3);

        // check for validity
        try {
            check_vector( v1 );
        } catch ( utility::excn::EXCN_Msg_Exception & e ) {
            std::string expected_error_message = "Unreasonable range for center or normal! Check your input vectors!";
            TS_ASSERT( expected_error_message == e.msg() );
        }
    }
    
    // average  embeddings
    void test_average_embeddings() {
        
        TS_TRACE("Test average embeddings");
        
        // define vectors
        Vector v1(1, 2, 3);
        Vector v2(4, 5, 6);
        Vector v3(7, 8, 9);
        Vector v4(7, 5, 3);
        Vector avg_center(4, 5, 6);
        Vector avg_normal(0.57735, 0.57735, 0.57735);
        
        // create embedding objects
        EmbeddingDefOP emb1( new EmbeddingDef( v1, v2 ) );
        EmbeddingDefOP emb2( new EmbeddingDef( v2, v3 ) );
        EmbeddingDefOP emb3( new EmbeddingDef( v3, v4 ) );
        
            // put them into a vector
        utility::vector1< EmbeddingDefOP > embeddings;
        embeddings.push_back( emb1 );
        embeddings.push_back( emb2 );
        embeddings.push_back( emb3 );
        
            // average them
        EmbeddingDefOP avg = average_embeddings( embeddings );
        
            // check for equality
        TS_ASSERT( position_equal_within_delta( avg->center(), avg_center, 0.001 ) );
        TS_ASSERT( position_equal_within_delta( avg->normal(), avg_normal, 0.001 ) );
    }

    // average antiparallel embeddings
    void test_average_antiparallel_embeddings() {
        
        TS_TRACE("Test average antiparallel embeddings");
        
        // define vectors
        Vector v1(1, 2, 3);
        Vector v2(4, 5, 6);
        Vector v3(1, 2, 9);
        Vector v4(3, 4, -8);
        Vector avg_center(2.5, 3.5, 4.5);
        Vector avg_normal(-0.116047, -0.116047, 0.986441);
        
        // create embedding objects
        EmbeddingDefOP emb1( new EmbeddingDef( v1, v3 ) );
        EmbeddingDefOP emb2( new EmbeddingDef( v2, v4 ) );
        
        // put them into a vector
        utility::vector1< EmbeddingDefOP > embeddings;
        embeddings.push_back( emb1 );
        embeddings.push_back( emb2 );
        
        // average them
        EmbeddingDefOP avg = average_antiparallel_embeddings( embeddings );
        
        // check for equality
        TS_ASSERT( position_equal_within_delta( avg->center(), avg_center, 0.001 ) );
        TS_ASSERT( position_equal_within_delta( avg->normal(), avg_normal, 0.001 ) );
    }
    
    // split topology by jump
    void test_split_topology_by_jump() {

        TS_TRACE("Test split topology by jump");
        using namespace core::conformation::membrane;
        using namespace protocols::membrane;
        using namespace protocols::membrane::geometry;

        // read in pose and create topology object
        Pose pose, pose_up, pose_down;
        core::import_pose::pose_from_pdb( pose, "protocols/membrane/geometry/1AFO_.pdb" );
        AddMembraneMoverOP addmem9( new AddMembraneMover( "protocols/membrane/geometry/1AFO__tr.span" ) );
        addmem9->apply(pose);
        
        SpanningTopology topo = SpanningTopology( "protocols/membrane/geometry/1AFO__tr.span" );
        SpanningTopology topo_up, topo_down;
        
        // call function
        split_topology_by_jump( pose, 1, topo, pose_up, pose_down, topo_up, topo_down );
        
        // test
        TS_ASSERT_EQUALS( topo_up.span(1)->start(), 15 );
        TS_ASSERT_EQUALS( topo_up.span(1)->end(), 31 );
        TS_ASSERT_EQUALS( topo_down.span(1)->start(), 15 );
        TS_ASSERT_EQUALS( topo_down.span(1)->end(), 33 );
    }

    // split topology by jump
    void test_split_topology_by_jump_noshift() {
        
        TS_TRACE("Test split topology by jump, no shift");
        using namespace core::conformation::membrane;
        using namespace protocols::membrane;
        using namespace protocols::membrane::geometry;
        
        // read in pose and create topology object
        Pose pose;
        core::import_pose::pose_from_pdb( pose, "protocols/membrane/geometry/1AFO_.pdb" );
        AddMembraneMoverOP addmem10( new AddMembraneMover( "protocols/membrane/geometry/1AFO__tr.span" ) );
        addmem10->apply(pose);
        
        SpanningTopologyOP topo( new SpanningTopology( "protocols/membrane/geometry/1AFO__tr.span" ) );
        SpanningTopologyOP topo_up( new SpanningTopology() );
        SpanningTopologyOP topo_down( new SpanningTopology() );
        
        // call function
        split_topology_by_jump_noshift( pose, 1, topo, topo_up, topo_down );
        
        // test
        TS_ASSERT_EQUALS( topo_up->span(1)->start(), 15 );
        TS_ASSERT_EQUALS( topo_up->span(1)->end(), 31 );
        TS_ASSERT_EQUALS( topo_down->span(1)->start(), 55 );
        TS_ASSERT_EQUALS( topo_down->span(1)->end(), 73 );
    }



private: 
	
	/// @brief Position equal within delta (helper method)
	bool position_equal_within_delta( Vector a, Vector b, Real delta ) {

		TS_ASSERT_DELTA( a.x(), b.x(), delta );
		TS_ASSERT_DELTA( a.y(), b.y(), delta );
		TS_ASSERT_DELTA( a.z(), b.z(), delta );

		return true;
	}

private:
	
	// Helix axes and foltrees
    PoseOP m2_pose_;
    PoseOP glpA_pose_;

	// RMSD calcs
    core::pose::PoseOP native_pose_;
    core::pose::PoseOP test_pose_;
    
}; // Protocols-level membrane utilities class
