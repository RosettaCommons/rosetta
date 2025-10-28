// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/membrane/MembraneUtil.cxxtest.hh
/// @brief   Unit Test: Protocols-level RosettaMP Framework calculations
///    Last Modified: 6/10/15
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

#include <core/kinematics/FoldTree.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/DockingPartners.hh>
#include <core/import_pose/import_pose.hh>
#include <core/types.hh>


// Utility Headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <string>

using namespace core;
using namespace core::import_pose;
using namespace core::conformation;
using namespace core::conformation::membrane;
using namespace core::pose;
using namespace core::kinematics;
using namespace protocols::membrane;
using namespace protocols::membrane::geometry;

static basic::Tracer TR("protocols.membrane.MembraneUtil.cxxtest");

/// @brief Unit Test suite for protocols-level membrane util class
class MembraneUtil : public CxxTest::TestSuite {

public: // test functions

	// Test Setup Functions ///////////////////////////

	/// @brief Setup the unit test
	void setUp() {

		//        using namespace core::import_pose;

		// Initialize core & options system
		core_init_with_additional_options( "-mp:lipids:has_pore 0" );
	}

	/// @brief Tear down the unit test
	void tearDown() {}

	// Test calc helix axis method
	void test_calc_helix_axis() {

		TR <<  "=========Testing the helix axes are correctly calculated for each span in the test set"  << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		// 1. TM domain of the M2 proton channel (single helix)
		PoseOP m2_pose ( new Pose() );
		pose_from_file( *m2_pose, "protocols/membrane/1mp6_transformed.pdb" , core::import_pose::PDB_file);
		AddMembraneMoverOP add_memb1 = utility::pointer::make_shared< AddMembraneMover >( "protocols/membrane/1mp6.span" );
		add_memb1->restore_lazaridis_IMM1_behavior( true );
		add_memb1->apply( *m2_pose );

		// 2. Glycophorin A (two helices, oriented 'somewhat' opposite one another)
		PoseOP glpA_pose( new Pose() );
		pose_from_file( *glpA_pose, "protocols/membrane/1AFO_AB.pdb" , core::import_pose::PDB_file);
		AddMembraneMoverOP add_memb2 = utility::pointer::make_shared< AddMembraneMover >( "protocols/membrane/1AFO_AB.span" );
		add_memb2->restore_lazaridis_IMM1_behavior( true );
		add_memb2->apply( *glpA_pose );

		// Calculate helix axis for each transembrane span
		Vector m2_first_span( calc_helix_axis( *m2_pose, 1 ) );
		Vector glpA_first_span( calc_helix_axis( *glpA_pose, 1 ) );
		Vector glpA_second_span( calc_helix_axis( *glpA_pose, 2 ) );

		Vector m2_first_expected( -3.208, 1.1903, 30.47033 );
		Vector glpA_first_expected( 1.409, -3.2603, 23.7866 );
		Vector glpA_second_expected( -0.886, 13.6927, 23.1553 );

		position_equal_within_delta( m2_first_span, m2_first_expected, 0.005 );
		position_equal_within_delta( glpA_first_span, glpA_first_expected, 0.005 );
		position_equal_within_delta( glpA_second_span, glpA_second_expected, 0.005 );

#endif

	}


	// Test Calc Angles Method

	/// @brief Calculate angles from single helix pose and two helix pose
	void test_calc_helix_angles() {

		TR <<  "=========Testing method for calculating the angle between a single helix and the membrane normal in a single helix pose"  << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		// 1. TM domain of the M2 proton channel (single helix)
		PoseOP m2_pose ( new Pose() );
		pose_from_file( *m2_pose, "protocols/membrane/1mp6_transformed.pdb" , core::import_pose::PDB_file);
		AddMembraneMoverOP add_memb1 = utility::pointer::make_shared< AddMembraneMover >( "protocols/membrane/1mp6.span" );
		add_memb1->restore_lazaridis_IMM1_behavior( true );
		add_memb1->apply( *m2_pose );

		// 2. Glycophorin A (two helices, oriented 'somewhat' opposite one another)
		PoseOP glpA_pose( new Pose() );
		pose_from_file( *glpA_pose, "protocols/membrane/1AFO_AB.pdb" , core::import_pose::PDB_file);
		AddMembraneMoverOP add_memb2 = utility::pointer::make_shared< AddMembraneMover >( "protocols/membrane/1AFO_AB.span" );
		add_memb2->restore_lazaridis_IMM1_behavior( true );
		add_memb2->apply( *glpA_pose );

		// Calculate helix tilt angles & compare
		Real m2_first_angle( calc_helix_tilt_angle( *m2_pose, 1 ) );
		Real glpA_first_angle( calc_helix_tilt_angle( *glpA_pose, 1 ) );
		Real glpA_second_angle( calc_helix_tilt_angle( *glpA_pose, 2 ) );

		TS_ASSERT_DELTA( m2_first_angle, 6.41, 0.005 );
		TS_ASSERT_DELTA( glpA_first_angle, 8.49, 0.005 );
		TS_ASSERT_DELTA( glpA_second_angle, 30.65, 0.005 );

#endif
	}


	// Test Angle RMSD method

	/// @brief Calculate RMSD between reference and new angle
	void test_angle_rmsd_method() {

		TR <<  "=========Testing method for calcuating the rms between a measured and reference angle value"  << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

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

#endif

	}

	///// These tests are separate and should be moved down into membrane info after debugging ////
	// There are more test cases here, but I'm only including the ones in use for now. Will extend later
	void test_is_fixed_on_moveable_memb() {

		TR <<  "=========Test is_membrane_fixed method on a moveable membrane"  << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		// 1. TM domain of the M2 proton channel (single helix)
		PoseOP m2_pose ( new Pose() );
		pose_from_file( *m2_pose, "protocols/membrane/1mp6_transformed.pdb" , core::import_pose::PDB_file);
		AddMembraneMoverOP add_memb1 = utility::pointer::make_shared< AddMembraneMover >( "protocols/membrane/1mp6.span" );
		add_memb1->apply( *m2_pose );

		// 2. Glycophorin A (two helices, oriented 'somewhat' opposite one another)
		PoseOP glpA_pose( new Pose() );
		pose_from_file( *glpA_pose, "protocols/membrane/1AFO_AB.pdb" , core::import_pose::PDB_file);
		AddMembraneMoverOP add_memb2 = utility::pointer::make_shared< AddMembraneMover >( "protocols/membrane/1AFO_AB.span" );
		add_memb2->apply( *glpA_pose );

		// Set up new foldtree for m2 rooted at the previous downstream residue
		// of the membrane jump
		Size m2_jump( m2_pose->conformation().membrane_info()->membrane_jump() );
		Size m2_downstream( m2_pose->conformation().fold_tree().downstream_jump_residue( m2_jump ) );
		FoldTreeOP m2_foldtree = utility::pointer::make_shared< FoldTree >( m2_pose->conformation().fold_tree() );
		m2_foldtree->reorder( m2_downstream );
		m2_pose->fold_tree( *m2_foldtree );
		m2_foldtree->show( TR );
		TS_ASSERT( !is_membrane_fixed( *m2_pose ) );

		// Set up new foldtree for glpA rooted at the previous downstream residue
		// of the membrane jump
		Size glpA_jump( glpA_pose->conformation().membrane_info()->membrane_jump() );
		Size glpA_downstream( glpA_pose->conformation().fold_tree().downstream_jump_residue( glpA_jump ) );
		FoldTreeOP glpA_foldtree = utility::pointer::make_shared< FoldTree >( glpA_pose->conformation().fold_tree() );
		glpA_foldtree->reorder( glpA_downstream );
		glpA_pose->fold_tree( *glpA_foldtree );
		glpA_foldtree->show( TR );
		TS_ASSERT( !is_membrane_fixed( *glpA_pose ) );

#endif

	}

	void test_is_fixed_on_fixed_memb() {

		TR <<  "=========Test is_fixed method on a fixed membrane"  << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		// 1. TM domain of the M2 proton channel (single helix)
		PoseOP m2_pose ( new Pose() );
		pose_from_file( *m2_pose, "protocols/membrane/1mp6_transformed.pdb" , core::import_pose::PDB_file);
		AddMembraneMoverOP add_memb1 = utility::pointer::make_shared< AddMembraneMover >( "protocols/membrane/1mp6.span" );
		add_memb1->apply( *m2_pose );

		// 2. Glycophorin A (two helices, oriented 'somewhat' opposite one another)
		PoseOP glpA_pose( new Pose() );
		pose_from_file( *glpA_pose, "protocols/membrane/1AFO_AB.pdb" , core::import_pose::PDB_file);
		AddMembraneMoverOP add_memb2 = utility::pointer::make_shared< AddMembraneMover >( "protocols/membrane/1AFO_AB.span" );
		add_memb2->apply( *glpA_pose );

		// When the pose is loaded in via AddMembraneMover, by default the foldtree
		// enforces a fixed membrane position
		TS_ASSERT( is_membrane_fixed( *m2_pose ) );
		TS_ASSERT( is_membrane_fixed( *glpA_pose ) );

#endif

	}

	void test_is_independently_moveable_on_fixed_memb() {

		TR <<  "=========Test is_independently_moveable on a fixed membrane"  << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		// 1. TM domain of the M2 proton channel (single helix)
		PoseOP m2_pose ( new Pose() );
		pose_from_file( *m2_pose, "protocols/membrane/1mp6_transformed.pdb" , core::import_pose::PDB_file);
		AddMembraneMoverOP add_memb1 = utility::pointer::make_shared< AddMembraneMover >( "protocols/membrane/1mp6.span" );
		add_memb1->apply( *m2_pose );

		// 2. Glycophorin A (two helices, oriented 'somewhat' opposite one another)
		PoseOP glpA_pose( new Pose() );
		pose_from_file( *glpA_pose, "protocols/membrane/1AFO_AB.pdb" , core::import_pose::PDB_file);
		AddMembraneMoverOP add_memb2 = utility::pointer::make_shared< AddMembraneMover >( "protocols/membrane/1AFO_AB.span" );
		add_memb2->apply( *glpA_pose );

		TS_ASSERT( !is_membrane_moveable_by_itself( *m2_pose ) );
		TS_ASSERT( !is_membrane_moveable_by_itself( *glpA_pose ) );

#endif

	}

	void test_is_independently_moveable_on_independently_moveable_memb() {

		TR <<  "=========Test is_independently_moveable on a moveable membrane"  << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		// 1. TM domain of the M2 proton channel (single helix)
		PoseOP m2_pose ( new Pose() );
		pose_from_file( *m2_pose, "protocols/membrane/1mp6_transformed.pdb" , core::import_pose::PDB_file);
		AddMembraneMoverOP add_memb1 = utility::pointer::make_shared< AddMembraneMover >( "protocols/membrane/1mp6.span" );
		add_memb1->apply( *m2_pose );

		// 2. Glycophorin A (two helices, oriented 'somewhat' opposite one another)
		PoseOP glpA_pose( new Pose() );
		pose_from_file( *glpA_pose, "protocols/membrane/1AFO_AB.pdb" , core::import_pose::PDB_file);
		AddMembraneMoverOP add_memb2 = utility::pointer::make_shared< AddMembraneMover >( "protocols/membrane/1AFO_AB.span" );
		add_memb2->apply( *glpA_pose );

		// Set up new foldtree for m2 rooted at the previous downstream residue
		// of the membrane jump
		Size m2_jump( m2_pose->conformation().membrane_info()->membrane_jump() );
		Size m2_downstream( m2_pose->conformation().fold_tree().downstream_jump_residue( m2_jump ) );
		FoldTreeOP m2_foldtree = utility::pointer::make_shared< FoldTree >( m2_pose->conformation().fold_tree() );
		m2_foldtree->reorder( m2_downstream );
		m2_pose->fold_tree( *m2_foldtree );
		m2_foldtree->show( TR );
		TS_ASSERT( is_membrane_moveable_by_itself( *m2_pose ) );

		// Set up new foldtree for glpA rooted at the previous downstream residue
		// of the membrane jump
		Size glpA_jump( glpA_pose->conformation().membrane_info()->membrane_jump() );
		Size glpA_downstream( glpA_pose->conformation().fold_tree().downstream_jump_residue( glpA_jump ) );
		FoldTreeOP glpA_foldtree = utility::pointer::make_shared< FoldTree >( glpA_pose->conformation().fold_tree() );
		glpA_foldtree->reorder( glpA_downstream );
		glpA_pose->fold_tree( *glpA_foldtree );
		glpA_foldtree->show( TR );
		TS_ASSERT( is_membrane_moveable_by_itself( *glpA_pose ) );

#endif

	}

	/// @brief Calculate membrane backbone rmsd with superposition
	void test_membrane_bb_rmsd_with_super() {

		TR <<  "=========Calculating membrane backbone rmsd with superposition"  << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		// Test Cases for membrane rmsd calculations
		// 3. Native Glycophprin A
		core::pose::PoseOP native_pose( new Pose() );
		pose_from_file( *native_pose, "protocols/membrane/1afo_in.pdb" , core::import_pose::PDB_file);

		// 4. Transformed glycophorin A
		core::pose::PoseOP test_pose( new Pose() );
		pose_from_file( *test_pose, "protocols/membrane/1afo_decoy.pdb" , core::import_pose::PDB_file);

		AddMembraneMoverOP add_memb3( new AddMembraneMover( "protocols/membrane/1afo_tr.span", 1 ) );
		add_memb3->apply( *native_pose );
		add_memb3->apply( *test_pose );

		core::Real rms = mem_bb_rmsd_with_super( *native_pose, *test_pose );
		TS_ASSERT_DELTA( rms, 0.0005, 0.0003 );


#endif
	}

	/// @brief Calculate membrane backbone rmsd with superposition
	void test_membrane_bb_rmsd_no_super() {

		TR <<  "=========Calculating membrane backbone rmsd without superposition"  << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		// Test Cases for membrane rmsd calculations
		// 3. Native Glycophprin A
		core::pose::PoseOP native_pose( new Pose() );
		pose_from_file( *native_pose, "protocols/membrane/1afo_in.pdb" , core::import_pose::PDB_file);

		// 4. Transformed glycophorin A
		core::pose::PoseOP test_pose( new Pose() );
		pose_from_file( *test_pose, "protocols/membrane/1afo_decoy.pdb" , core::import_pose::PDB_file);

		AddMembraneMoverOP add_memb3( new AddMembraneMover( "protocols/membrane/1afo_tr.span", 1 ) );
		add_memb3->apply( *native_pose );
		add_memb3->apply( *test_pose );

		core::Real rms = mem_bb_rmsd_no_super( *native_pose, *test_pose );
		TS_ASSERT_DELTA( rms, 3.8042, 0.0003 );

#endif

	}

	/// @brief Calculate membrane all atom rmsd with superposition
	void test_membrane_bb_rmsd_with_super_allatom() {

		TR <<  "=========Calculating membrane allatom rmsd with superposition"  << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		// Test Cases for membrane rmsd calculations
		// 3. Native Glycophprin A
		core::pose::PoseOP native_pose( new Pose() );
		pose_from_file( *native_pose, "protocols/membrane/1afo_in.pdb" , core::import_pose::PDB_file);

		// 4. Transformed glycophorin A
		core::pose::PoseOP test_pose( new Pose() );
		pose_from_file( *test_pose, "protocols/membrane/1afo_decoy.pdb" , core::import_pose::PDB_file);

		AddMembraneMoverOP add_memb3( new AddMembraneMover( "protocols/membrane/1afo_tr.span", 1 ) );
		add_memb3->apply( *native_pose );
		add_memb3->apply( *test_pose );

		core::Real rms = mem_all_atom_rmsd_with_super( *native_pose, *test_pose );
		TS_ASSERT_DELTA( rms, 0.5052, 0.0003 );

#endif

	}

	/// @brief Calculate membrane all atom rmsd without superposition
	void test_membrane_bb_rmsd_no_super_allatom() {

		TR <<  "=========Calculating membrane backbone rmsd without superposition"  << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		// Test Cases for membrane rmsd calculations
		// 3. Native Glycophprin A
		core::pose::PoseOP native_pose( new Pose() );
		pose_from_file( *native_pose, "protocols/membrane/1afo_in.pdb" , core::import_pose::PDB_file);

		// 4. Transformed glycophorin A
		core::pose::PoseOP test_pose( new Pose() );
		pose_from_file( *test_pose, "protocols/membrane/1afo_decoy.pdb" , core::import_pose::PDB_file);

		AddMembraneMoverOP add_memb3( new AddMembraneMover( "protocols/membrane/1afo_tr.span", 1 ) );
		add_memb3->apply( *native_pose );
		add_memb3->apply( *test_pose );

		core::Real rms = mem_all_atom_rmsd_no_super( *native_pose, *test_pose );
		TS_ASSERT_DELTA( rms, 3.9376, 0.0003 );

#endif

	}

	// compute_structure_based_embedding
	void test_compute_structure_based_embedding1() {

		TR << "=========Test compute_structure_based_embedding 1" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION


		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;

		// 1AFO
		// read in pose and create topology object
		TR << "1AFO" << std::endl;
		Pose pose1;
		core::import_pose::pose_from_file( pose1, "protocols/membrane/geometry/1AFO_.pdb" , core::import_pose::PDB_file);

		// create membrane pose
		AddMembraneMoverOP addmem1( new AddMembraneMover( "protocols/membrane/geometry/1AFO__tr.span" ) );
		addmem1->restore_lazaridis_IMM1_behavior( true );
		addmem1->apply( pose1 );

		// define vectors and object
		Vector center1(0.53225, 0.361, 0.095);
		Vector normal1(-0.849116, -0.512025, 0.129741);

		// compute embedding
		EmbeddingDefOP embed1( compute_structure_based_embedding( pose1 ) );

		// compare
		TS_ASSERT( position_equal_within_delta( embed1->center(), center1, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( embed1->normal(), normal1, 0.001 ) );

#endif

	}

	// compute_structure_based_embedding
	void test_compute_structure_based_embedding2() {

		TR << "=========Test compute_structure_based_embedding 2" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;

		// 1BL8
		TR << "1BL8" << std::endl;
		Vector center( 0, 0, 0 );
		Vector normal( 0, 0, 1 );
		Pose pose2;
		core::import_pose::pose_from_file( pose2, "protocols/membrane/geometry/1BL8_.pdb" , core::import_pose::PDB_file);
		AddMembraneMoverOP addmem2( new AddMembraneMover( "protocols/membrane/geometry/1BL8__tr.span" ) );
		addmem2->restore_lazaridis_IMM1_behavior( true );
		addmem2->apply(pose2);
		Vector center2(73.9421, 26.7549, 24.4493);
		Vector normal2(0.384026, -0.0403822, 0.922439);
		EmbeddingDefOP embed2( compute_structure_based_embedding( pose2 ) );
		TS_ASSERT( position_equal_within_delta( embed2->center(), center2, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( embed2->normal(), normal2, 0.001 ) );

#endif

	}

	// compute_structure_based_embedding
	void test_compute_structure_based_embedding3() {

		TR << "=========Test compute_structure_based_embedding 3" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;

		// 1QJP - beta-barrel
		TR << "1QJP" << std::endl;
		Vector center( 0, 0, 0 );
		Vector normal( 0, 0, 1 );
		Pose pose3;
		core::import_pose::pose_from_file( pose3, "protocols/membrane/geometry/1QJP_.pdb" , core::import_pose::PDB_file);
		AddMembraneMoverOP addmem3( new AddMembraneMover( "protocols/membrane/geometry/1QJP__tr.span" ) );
		addmem3->restore_lazaridis_IMM1_behavior( true );
		addmem3->apply(pose3);
		Vector center3(31.2161, 16.9685, 37.6119);
		Vector normal3(0.877926, -0.47167, 0.0822946);
		EmbeddingDefOP embed3( compute_structure_based_embedding( pose3 ) );
		TS_ASSERT( position_equal_within_delta( embed3->center(), center3, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( embed3->normal(), normal3, 0.001 ) );

#endif

	}

	// compute_structure_based_embedding
	void test_compute_structure_based_embedding4() {

		TR << "=========Test compute_structure_based_embedding 4" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;

		// 2BS2
		TR << "2BS2" << std::endl;
		Vector center( 0, 0, 0 );
		Vector normal( 0, 0, 1 );
		Pose pose4;
		core::import_pose::pose_from_file( pose4, "protocols/membrane/geometry/2BS2_CF.pdb" , core::import_pose::PDB_file);
		AddMembraneMoverOP addmem4( new AddMembraneMover( "protocols/membrane/geometry/2BS2_CF_tr.span" ) );
		addmem4->restore_lazaridis_IMM1_behavior( true );
		addmem4->apply(pose4);
		Vector center4(21.4326, 6.0464, -41.0573);
		Vector normal4(0.0060, 0.0117348, 0.9999);
		EmbeddingDefOP embed4( compute_structure_based_embedding( pose4 ) );
		TS_ASSERT( position_equal_within_delta( embed4->center(), center4, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( embed4->normal(), normal4, 0.001 ) );

#endif

	}

	// compute_structure_based_embedding
	void test_compute_structure_based_embedding5() {

		TR << "=========Test compute_structure_based_embedding 5" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;

		// 2MPN
		TR << "2MPN" << std::endl;
		Vector center( 0, 0, 0 );
		Vector normal( 0, 0, 1 );
		Pose pose5;
		core::import_pose::pose_from_file( pose5, "protocols/membrane/geometry/2MPN_.pdb" , core::import_pose::PDB_file);
		AddMembraneMoverOP addmem5( new AddMembraneMover( "protocols/membrane/geometry/2MPN__tr.span" ) );
		addmem5->restore_lazaridis_IMM1_behavior( true );
		addmem5->apply(pose5);
		Vector center5(0.3645, 3.66025, 41.345);
		Vector normal5(0.006956, 0.996657, 0.0813998);
		EmbeddingDefOP embed5( compute_structure_based_embedding( pose5 ) );
		TS_ASSERT( position_equal_within_delta( embed5->center(), center5, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( embed5->normal(), normal5, 0.001 ) );

#endif

	}

	// compute_structure_based_embedding
	void test_compute_structure_based_embedding6() {

		TR << "=========Test compute_structure_based_embedding 6" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;

		// 2OAR
		TR << "2OAR" << std::endl;
		Vector center( 0, 0, 0 );
		Vector normal( 0, 0, 1 );
		Pose pose6;
		core::import_pose::pose_from_file( pose6, "protocols/membrane/geometry/2OAR_.pdb" , core::import_pose::PDB_file);
		AddMembraneMoverOP addmem6( new AddMembraneMover( "protocols/membrane/geometry/2OAR__tr.span" ) );
		addmem6->restore_lazaridis_IMM1_behavior( true );
		addmem6->apply(pose6);
		Vector center6(18.8453, 122.117, 1.079);
		Vector normal6(-0.755094, -0.655578, 0.00710981);
		EmbeddingDefOP embed6( compute_structure_based_embedding( pose6 ) );
		TS_ASSERT( position_equal_within_delta( embed6->center(), center6, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( embed6->normal(), normal6, 0.001 ) );

#endif

	}

	// compute_structure_based_embedding
	void test_compute_structure_based_embedding7() {

		TR << "=========Test compute_structure_based_embedding 7" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION


		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;

		// 2UUH
		TR << "2UUH" << std::endl;
		Vector center( 0, 0, 0 );
		Vector normal( 0, 0, 1 );
		Pose pose7;
		core::import_pose::pose_from_file( pose7, "protocols/membrane/geometry/2UUH__tr.pdb" , core::import_pose::PDB_file);
		AddMembraneMoverOP addmem7( new AddMembraneMover( "protocols/membrane/geometry/2UUH__tr.span" ) );
		addmem7->restore_lazaridis_IMM1_behavior( true );
		addmem7->apply(pose7);
		Vector center7(-0.000166667, -0.000125, 0.295625);
		Vector normal7(0, 0, 1);
		EmbeddingDefOP embed7( compute_structure_based_embedding( pose7 ) );
		TS_ASSERT( position_equal_within_delta( embed7->center(), center7, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( embed7->normal(), normal7, 0.001 ) );

#endif

	}

	// compute_structure_based_embedding
	void test_compute_structure_based_embedding8() {

		TR << "=========Test compute_structure_based_embedding 8" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;

		// 3PXO
		TR << "3PXO" << std::endl;
		Vector center( 0, 0, 0 );
		Vector normal( 0, 0, 1 );
		Pose pose8;
		core::import_pose::pose_from_file( pose8, "protocols/membrane/geometry/3PXO_.pdb" , core::import_pose::PDB_file);
		AddMembraneMoverOP addmem8( new AddMembraneMover( "protocols/membrane/geometry/3PXO__tr.span" ) );
		addmem8->restore_lazaridis_IMM1_behavior( true );
		addmem8->apply(pose8);
		Vector center8(-36.1201, -7.59636, 37.6713);
		Vector normal8(-0.9862, -0.164797, 0.0158378);
		EmbeddingDefOP embed8( compute_structure_based_embedding( pose8 ) );
		TS_ASSERT( position_equal_within_delta( embed8->center(), center8, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( embed8->normal(), normal8, 0.001 ) );

#endif

	}

	// check vector for reasonable size
	void test_check_vector() {

		TR << "=========Test check vector" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		// define vectors and object
		Vector v1(1, 2, 1000);
		Vector v2(4, 5000, 3);

		// check for validity
		try {
			check_vector( v1 );
		} catch (utility::excn::Exception & e ) {
			std::string expected_error_message = "Unreasonable range for center or normal! Check your input vectors!";
			TS_ASSERT( e.msg().find(expected_error_message) != std::string::npos );
		}

#endif
	}

	// average  embeddings
	void test_average_embeddings() {

		TR << "=========Test average embeddings" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

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

#endif

	}

	// average antiparallel embeddings
	void test_average_antiparallel_embeddings() {

		TR << "=========Test average antiparallel embeddings" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

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

#endif

	}

	// split topology by jump
	void test_split_topology_by_jump() {

		TR << "=========Test split topology by jump" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		// read in pose and create topology object
		Pose pose, pose_up, pose_down;
		core::import_pose::pose_from_file( pose, "protocols/membrane/geometry/1AFO_.pdb" , core::import_pose::PDB_file);
		AddMembraneMoverOP addmem9( new AddMembraneMover( "protocols/membrane/geometry/1AFO__tr.span" ) );
		addmem9->apply(pose);

		std::map< std::string, core::Size > pdb2pose_map = core::pose::get_pdb2pose_numbering_as_stdmap( pose );
		SpanningTopology topo = SpanningTopology( "protocols/membrane/geometry/1AFO__tr.span", pdb2pose_map );
		SpanningTopology topo_up, topo_down;

		// call function
		split_topology_by_jump( pose, 1, topo, pose_up, pose_down, topo_up, topo_down );

		// test
		TS_ASSERT_EQUALS( topo_up.span(1)->start(), 15 );
		TS_ASSERT_EQUALS( topo_up.span(1)->end(), 31 );
		TS_ASSERT_EQUALS( topo_down.span(1)->start(), 15 );
		TS_ASSERT_EQUALS( topo_down.span(1)->end(), 33 );

#endif
	}

	// split topology by jump
	void test_split_topology_by_jump_noshift() {

		TR << "=========Test split topology by jump, no shift" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		// read in pose and create topology object
		Pose pose;
		core::import_pose::pose_from_file( pose, "protocols/membrane/geometry/1AFO_.pdb" , core::import_pose::PDB_file);
		AddMembraneMoverOP addmem10( new AddMembraneMover( "protocols/membrane/geometry/1AFO__tr.span" ) );
		addmem10->apply(pose);

		std::map< std::string, core::Size > pdb2pose_map = core::pose::get_pdb2pose_numbering_as_stdmap( pose );
		SpanningTopologyOP topo( new SpanningTopology( "protocols/membrane/geometry/1AFO__tr.span", pdb2pose_map ) );
		SpanningTopologyOP topo_up( new SpanningTopology() );
		SpanningTopologyOP topo_down( new SpanningTopology() );

		// call function
		split_topology_by_jump_noshift( pose, 1, topo, topo_up, topo_down );

		// test
		TS_ASSERT_EQUALS( topo_up->span(1)->start(), 15 );
		TS_ASSERT_EQUALS( topo_up->span(1)->end(), 31 );
		TS_ASSERT_EQUALS( topo_down->span(1)->start(), 55 );
		TS_ASSERT_EQUALS( topo_down->span(1)->end(), 73 );

#endif

	}

	// split topology by chain
	void test_split_topology_by_chain_noshift() {

		TR << "=========Test split topology by chain" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		// read in pose and create topology object
		Pose pose;
		core::import_pose::pose_from_file( pose, "protocols/membrane/1C17_tr.pdb" , core::import_pose::PDB_file);
		std::map< std::string, core::Size > pdb2pose_map = core::pose::get_pdb2pose_numbering_as_stdmap( pose );
		SpanningTopologyOP topo( new SpanningTopology( "protocols/membrane/1C17_tr.span", pdb2pose_map ) );

		// call function
		utility::vector1< SpanningTopologyOP > topos( split_topology_by_chain_noshift( pose, topo ) );

		// test
		TS_ASSERT_EQUALS( topos[1]->span(1)->start(), 12 );
		TS_ASSERT_EQUALS( topos[1]->span(1)->end(), 32 );
		TS_ASSERT_EQUALS( topos[1]->span(2)->start(), 52 );
		TS_ASSERT_EQUALS( topos[1]->span(2)->end(), 71 );

		TS_ASSERT_EQUALS( topos[2]->span(1)->start(), 91 );
		TS_ASSERT_EQUALS( topos[2]->span(1)->end(), 111 );
		TS_ASSERT_EQUALS( topos[2]->span(2)->start(), 131 );
		TS_ASSERT_EQUALS( topos[2]->span(2)->end(), 150 );

		TS_ASSERT_EQUALS( topos[3]->span(1)->start(), 170 );
		TS_ASSERT_EQUALS( topos[3]->span(1)->end(), 190 );
		TS_ASSERT_EQUALS( topos[3]->span(2)->start(), 210 );
		TS_ASSERT_EQUALS( topos[3]->span(2)->end(), 229 );

		TS_ASSERT_EQUALS( topos[4]->span(1)->start(), 249 );
		TS_ASSERT_EQUALS( topos[4]->span(1)->end(), 269 );
		TS_ASSERT_EQUALS( topos[4]->span(2)->start(), 289 );
		TS_ASSERT_EQUALS( topos[4]->span(2)->end(), 308 );

		TS_ASSERT_EQUALS( topos[5]->span(1)->start(), 328 );
		TS_ASSERT_EQUALS( topos[5]->span(1)->end(), 348 );
		TS_ASSERT_EQUALS( topos[5]->span(2)->start(), 368 );
		TS_ASSERT_EQUALS( topos[5]->span(2)->end(), 387 );

		TS_ASSERT_EQUALS( topos[6]->span(1)->start(), 407 );
		TS_ASSERT_EQUALS( topos[6]->span(1)->end(), 427 );
		TS_ASSERT_EQUALS( topos[6]->span(2)->start(), 447 );
		TS_ASSERT_EQUALS( topos[6]->span(2)->end(), 466 );

		TS_ASSERT_EQUALS( topos[7]->span(1)->start(), 486 );
		TS_ASSERT_EQUALS( topos[7]->span(1)->end(), 506 );
		TS_ASSERT_EQUALS( topos[7]->span(2)->start(), 526 );
		TS_ASSERT_EQUALS( topos[7]->span(2)->end(), 545 );

		TS_ASSERT_EQUALS( topos[8]->span(1)->start(), 565 );
		TS_ASSERT_EQUALS( topos[8]->span(1)->end(), 585 );
		TS_ASSERT_EQUALS( topos[8]->span(2)->start(), 605 );
		TS_ASSERT_EQUALS( topos[8]->span(2)->end(), 624 );

		TS_ASSERT_EQUALS( topos[9]->span(1)->start(), 644 );
		TS_ASSERT_EQUALS( topos[9]->span(1)->end(), 664 );
		TS_ASSERT_EQUALS( topos[9]->span(2)->start(), 684 );
		TS_ASSERT_EQUALS( topos[9]->span(2)->end(), 703 );

		TS_ASSERT_EQUALS( topos[10]->span(1)->start(), 723 );
		TS_ASSERT_EQUALS( topos[10]->span(1)->end(), 744 );
		TS_ASSERT_EQUALS( topos[10]->span(2)->start(), 763 );
		TS_ASSERT_EQUALS( topos[10]->span(2)->end(), 782 );

		TS_ASSERT_EQUALS( topos[11]->span(1)->start(), 802 );
		TS_ASSERT_EQUALS( topos[11]->span(1)->end(), 823 );
		TS_ASSERT_EQUALS( topos[11]->span(2)->start(), 842 );
		TS_ASSERT_EQUALS( topos[11]->span(2)->end(), 861 );

		TS_ASSERT_EQUALS( topos[12]->span(1)->start(), 881 );
		TS_ASSERT_EQUALS( topos[12]->span(1)->end(), 901 );
		TS_ASSERT_EQUALS( topos[12]->span(2)->start(), 920 );
		TS_ASSERT_EQUALS( topos[12]->span(2)->end(), 940 );

		TS_ASSERT_EQUALS( topos[13]->span(1)->start(), 956 );
		TS_ASSERT_EQUALS( topos[13]->span(1)->end(), 975 );
		TS_ASSERT_EQUALS( topos[13]->span(2)->start(), 999 );
		TS_ASSERT_EQUALS( topos[13]->span(2)->end(), 1017 );
		TS_ASSERT_EQUALS( topos[13]->span(3)->start(), 1027 );
		TS_ASSERT_EQUALS( topos[13]->span(3)->end(), 1044 );
		TS_ASSERT_EQUALS( topos[13]->span(4)->start(), 1068 );
		TS_ASSERT_EQUALS( topos[13]->span(4)->end(), 1087 );

#endif

	}

	// compute embeddings by chain
	void test_compute_embeddings_by_chain() {

		TR << "=========Test compute embeddings by chain" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		// read in pose and create topology object
		Pose pose;
		core::import_pose::pose_from_file( pose, "protocols/membrane/1C17_tr.pdb" , core::import_pose::PDB_file);

		// call AddMembraneMover for spanfile
		AddMembraneMoverOP addmem( new AddMembraneMover( "protocols/membrane/1C17_tr.span" ) );
		addmem->apply(pose);

		// compute embeddings by chain
		EmbeddingOP embeddings( compute_embeddings_by_chain( pose ) );

		// initialize vectors
		core::Vector center1( 18.882, -0.12625, -0.28525 );
		core::Vector normal1( 0.064413, 0.0799618, 0.994715 );
		core::Vector center2( 18.6088, 11.9908, -0.248 );
		core::Vector normal2( 0.00767857, 0.0858815, 0.996276 );
		core::Vector center3( 12.5838, 22.4555, -0.05 );
		core::Vector normal3( -0.0432722, 0.0666127, 0.99684 );
		core::Vector center4( 2.26575, 28.673, 0.13175 );
		core::Vector normal4( -0.0763435, 0.0286892, 0.996669 );
		core::Vector center5( -9.813, 28.848, 0.24375 );
		core::Vector normal5( -0.0821136, -0.0251727, 0.996305 );
		core::Vector center6( -20.1778, 22.7, 0.26525 );
		core::Vector normal6( -0.0632901, -0.0720583, 0.99539 );
		core::Vector center7( -25.7978, 11.8975, 0.219 );
		core::Vector normal7( -0.0292156, -0.0935649, 0.995184 );
		core::Vector center8( -25.2455, -0.2675, 0.196 );
		core::Vector normal8( 0.0121257, -0.0949157, 0.995411 );
		core::Vector center9( -18.77, -10.611, 0.31275 );
		core::Vector normal9( 0.0382791, -0.0768116, 0.996311 );
		core::Vector center10( -8.5835, -16.2085, 0.07625 );
		core::Vector normal10( 0.0677452, -0.09455, 0.993212 );
		core::Vector center11( 3.07125, -17.2972, 0.29 );
		core::Vector normal11( 0.105914, -0.0375546, 0.993666 );
		core::Vector center12( 16.2923, -13.9145, 0.3135 );
		core::Vector normal12( -0.0246791, 0.175757, 0.984124 );
		core::Vector center13( 17.5571, -34.6114, 0.1045 );
		core::Vector normal13( 0.0546472, 0.0975181, 0.993732 );

		// test
		position_equal_within_delta( embeddings->embedding(1)->center(), center1 , 0.001 );
		position_equal_within_delta( embeddings->embedding(1)->normal(), normal1 , 0.001 );
		position_equal_within_delta( embeddings->embedding(2)->center(), center2 , 0.001 );
		position_equal_within_delta( embeddings->embedding(2)->normal(), normal2 , 0.001 );
		position_equal_within_delta( embeddings->embedding(3)->center(), center3 , 0.001 );
		position_equal_within_delta( embeddings->embedding(3)->normal(), normal3 , 0.001 );
		position_equal_within_delta( embeddings->embedding(4)->center(), center4 , 0.001 );
		position_equal_within_delta( embeddings->embedding(4)->normal(), normal4 , 0.001 );
		position_equal_within_delta( embeddings->embedding(5)->center(), center5 , 0.001 );
		position_equal_within_delta( embeddings->embedding(5)->normal(), normal5 , 0.001 );
		position_equal_within_delta( embeddings->embedding(6)->center(), center6 , 0.001 );
		position_equal_within_delta( embeddings->embedding(6)->normal(), normal6 , 0.001 );
		position_equal_within_delta( embeddings->embedding(7)->center(), center7 , 0.001 );
		position_equal_within_delta( embeddings->embedding(7)->normal(), normal7 , 0.001 );
		position_equal_within_delta( embeddings->embedding(8)->center(), center8 , 0.001 );
		position_equal_within_delta( embeddings->embedding(8)->normal(), normal8 , 0.001 );
		position_equal_within_delta( embeddings->embedding(9)->center(), center9 , 0.001 );
		position_equal_within_delta( embeddings->embedding(9)->normal(), normal9 , 0.001 );
		position_equal_within_delta( embeddings->embedding(10)->center(), center10 , 0.001 );
		position_equal_within_delta( embeddings->embedding(10)->normal(), normal10 , 0.001 );
		position_equal_within_delta( embeddings->embedding(11)->center(), center11 , 0.001 );
		position_equal_within_delta( embeddings->embedding(11)->normal(), normal11 , 0.001 );
		position_equal_within_delta( embeddings->embedding(12)->center(), center12 , 0.001 );
		position_equal_within_delta( embeddings->embedding(12)->normal(), normal12 , 0.001 );
		position_equal_within_delta( embeddings->embedding(13)->center(), center13 , 0.001 );
		position_equal_within_delta( embeddings->embedding(13)->normal(), normal13 , 0.001 );

#endif

	}

	// chain center-of-mass
	void test_chain_com() {

		TR << "=========Test chain center-of-mass" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		using namespace core::conformation::membrane;
		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;

		// read in pose and create topology object
		Pose pose;
		core::import_pose::pose_from_file( pose, "protocols/membrane/3EFF_TR.pdb" , core::import_pose::PDB_file);

		// call AddMembraneMover for spanfile
		AddMembraneMoverOP addmem( new AddMembraneMover( "protocols/membrane/3EFF_TR.span" ) );
		addmem->apply(pose);

		// get chains
		utility::vector1< int > chains = get_chains( pose );

		// initializations; hpark: slight chance with beta_nov15
		core::Vector com1( /*-8.6927, 3.7495, 19.4293   */-8.66951, 3.74931, 19.4430 );
		core::Vector com2( /*-4.70189, -9.09737, 19.3467*/-4.70721, -9.09404, 19.3549);
		core::Vector com3( /*8.92337, -3.8958, 18.6790  */ 8.90770,-3.89594, 18.7045 );
		core::Vector com4( /*4.4793, 9.22617, 19.3396   */ 4.46862, 9.22600, 19.3294 );

		// compute chain COM for all chains and pushback into vector
		core::Vector com_1 = protocols::membrane::chain_com( pose, chains[ 1 ] );
		core::Vector com_2 = protocols::membrane::chain_com( pose, chains[ 2 ] );
		core::Vector com_3 = protocols::membrane::chain_com( pose, chains[ 3 ] );
		core::Vector com_4 = protocols::membrane::chain_com( pose, chains[ 4 ] );

		TR << "core::Vector com1( " << com_1[0] << ", " << com_1[1] << ", " << com_1[2] << " );" << std::endl;
		TR << "core::Vector com2( " << com_2[0] << ", " << com_2[1] << ", " << com_2[2] << " );" << std::endl;
		TR << "core::Vector com3( " << com_3[0] << ", " << com_3[1] << ", " << com_3[2] << " );" << std::endl;
		TR << "core::Vector com4( " << com_4[0] << ", " << com_4[1] << ", " << com_4[2] << " );" << std::endl;

		// test
		position_equal_within_delta( com1, com_1 , 0.001 );
		position_equal_within_delta( com2, com_2 , 0.001 );
		position_equal_within_delta( com3, com_3 , 0.001 );
		position_equal_within_delta( com4, com_4 , 0.001 );

#endif

	}

	// residue closest to chain center-of-mass
	void test_rsd_closest_chain_com() {

		TR << "=========Test residue closest to chain center-of-mass" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		using namespace core::conformation::membrane;
		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;

		// read in pose and create topology object
		Pose pose;
		core::import_pose::pose_from_file( pose, "protocols/membrane/3EFF_TR.pdb" , core::import_pose::PDB_file);

		// call AddMembraneMover for spanfile
		AddMembraneMoverOP addmem( new AddMembraneMover( "protocols/membrane/3EFF_TR.span" ) );
		addmem->apply(pose);

		// get chains
		utility::vector1< int > chains = get_chains( pose );

		// get residue closest to COM
		core::Size closest1 = rsd_closest_to_chain_com( pose, chains[ 1 ] );
		core::Size closest2 = rsd_closest_to_chain_com( pose, chains[ 2 ] );
		core::Size closest3 = rsd_closest_to_chain_com( pose, chains[ 3 ] );
		core::Size closest4 = rsd_closest_to_chain_com( pose, chains[ 4 ] );

		// test
		TS_ASSERT_EQUALS( closest1, 88 );
		TS_ASSERT_EQUALS( closest2, 227 );
		TS_ASSERT_EQUALS( closest3, 366 );
		TS_ASSERT_EQUALS( closest4, 505 );

#endif
	}

	// per chain TM center-of-mass
	void test_tm_com() {

		TR << "=========Test per chain transmembrane center-of-mass" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		using namespace core::conformation::membrane;
		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;

		// read in pose and create topology object
		Pose pose;
		core::import_pose::pose_from_file( pose, "protocols/membrane/3EFF_TR.pdb" , core::import_pose::PDB_file);

		// call AddMembraneMover for spanfile
		AddMembraneMoverOP addmem( new AddMembraneMover( "protocols/membrane/3EFF_TR.span" ) );
		addmem->apply(pose);

		// initializations
		core::Vector tm1center( -10.4263, 9.96149, -2.58642 );
		core::Vector tm2center( -11.0477, -9.70225, -2.83403 );
		core::Vector tm3center( 9.87054, -9.53731, -2.45051 );
		core::Vector tm4center( 10.2604, 10.7574, -2.62829 );

		// get chains
		utility::vector1< int > chains = get_chains( pose );

		// compute COMs for TMspans for each chain
		core::Vector tmcom1 = chain_tm_com( pose, chains[ 1 ] );
		core::Vector tmcom2 = chain_tm_com( pose, chains[ 2 ] );
		core::Vector tmcom3 = chain_tm_com( pose, chains[ 3 ] );
		core::Vector tmcom4 = chain_tm_com( pose, chains[ 4 ] );

		// test
		position_equal_within_delta( tmcom1, tm1center , 0.001 );
		position_equal_within_delta( tmcom2, tm2center , 0.001 );
		position_equal_within_delta( tmcom3, tm3center , 0.001 );
		position_equal_within_delta( tmcom4, tm4center , 0.001 );

#endif
	}

	// residue closest to chain tm center-of-mass
	void test_rsd_closest_chain_tm_com() {

		TR << "=========Test residue closest to chain TM center-of-mass" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		using namespace core::conformation::membrane;
		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;

		// read in pose and create topology object
		Pose pose;
		core::import_pose::pose_from_file( pose, "protocols/membrane/3EFF_TR.pdb" , core::import_pose::PDB_file);

		// call AddMembraneMover for spanfile
		AddMembraneMoverOP addmem( new AddMembraneMover( "protocols/membrane/3EFF_TR.span" ) );
		addmem->apply(pose);

		// get chains
		utility::vector1< int > chains = get_chains( pose );

		// get residue closest to COM
		core::Size closest1 = rsd_closest_to_chain_tm_com( pose, chains[ 1 ] );
		core::Size closest2 = rsd_closest_to_chain_tm_com( pose, chains[ 2 ] );
		core::Size closest3 = rsd_closest_to_chain_tm_com( pose, chains[ 3 ] );
		core::Size closest4 = rsd_closest_to_chain_tm_com( pose, chains[ 4 ] );

		// test
		TS_ASSERT_EQUALS( closest1, 19 );
		TS_ASSERT_EQUALS( closest2, 158 );
		TS_ASSERT_EQUALS( closest3, 352 );
		TS_ASSERT_EQUALS( closest4, 491 );

#endif

	} // rsd closest to chain TM COM

	// create membrane foldtree anchor COM
	void test_create_membrane_foldtree_anchor_com() {

		TR << "=========Test create membrane foldtree anchor at center-of-mass" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		using namespace core::conformation::membrane;
		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;

		// read in pose and create topology object
		Pose pose;
		core::import_pose::pose_from_file( pose, "protocols/membrane/3EFF_TR.pdb" , core::import_pose::PDB_file);

		// call AddMembraneMover for spanfile
		AddMembraneMoverOP addmem( new AddMembraneMover( "protocols/membrane/3EFF_TR.span" ) );
		addmem->apply(pose);

		// create foldtree
		create_membrane_foldtree_anchor_com( pose );

		// test
		TS_ASSERT_EQUALS( pose.fold_tree().upstream_jump_residue( 1 ), 557 );
		TS_ASSERT_EQUALS( pose.fold_tree().downstream_jump_residue( 1 ), 88 );
		TS_ASSERT_EQUALS( pose.fold_tree().upstream_jump_residue( 2 ), 88 );
		TS_ASSERT_EQUALS( pose.fold_tree().downstream_jump_residue( 2 ), 227 );
		TS_ASSERT_EQUALS( pose.fold_tree().upstream_jump_residue( 3 ), 88 );
		TS_ASSERT_EQUALS( pose.fold_tree().downstream_jump_residue( 3 ), 366 );
		TS_ASSERT_EQUALS( pose.fold_tree().upstream_jump_residue( 4 ), 88 );
		TS_ASSERT_EQUALS( pose.fold_tree().downstream_jump_residue( 4 ), 505 );

#endif

	} // create membrane foldtree anchor COM

	// create membrane foldtree anchor TM COM
	void test_create_membrane_foldtree_anchor_tmcom() {

		TR << "=========Test create membrane foldtree anchor at TM center-of-mass" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		using namespace core::conformation::membrane;
		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;

		// read in pose and create topology object
		Pose pose;
		core::import_pose::pose_from_file( pose, "protocols/membrane/3EFF_TR.pdb" , core::import_pose::PDB_file);

		// call AddMembraneMover for spanfile
		AddMembraneMoverOP addmem( new AddMembraneMover( "protocols/membrane/3EFF_TR.span" ) );
		addmem->apply(pose);

		// create foldtree
		create_membrane_foldtree_anchor_tmcom( pose );

		// test
		TS_ASSERT_EQUALS( pose.fold_tree().upstream_jump_residue( 1 ), 557 );
		TS_ASSERT_EQUALS( pose.fold_tree().downstream_jump_residue( 1 ), 19 );
		TS_ASSERT_EQUALS( pose.fold_tree().upstream_jump_residue( 2 ), 19 );
		TS_ASSERT_EQUALS( pose.fold_tree().downstream_jump_residue( 2 ), 158 );
		TS_ASSERT_EQUALS( pose.fold_tree().upstream_jump_residue( 3 ), 19 );
		TS_ASSERT_EQUALS( pose.fold_tree().downstream_jump_residue( 3 ), 352 );
		TS_ASSERT_EQUALS( pose.fold_tree().upstream_jump_residue( 4 ), 19 );
		TS_ASSERT_EQUALS( pose.fold_tree().downstream_jump_residue( 4 ), 491 );

#endif
	} // create membrane foldtree anchor TM COM

	// pose TM COM
	void test_pose_tm_COM() {

		TR << "=========Test pose TM center-of-mass" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		using namespace core::conformation::membrane;
		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;

		// read in pose and create topology object
		Pose pose;
		core::import_pose::pose_from_file( pose, "protocols/membrane/3EFF_TR.pdb" , core::import_pose::PDB_file);

		// call AddMembraneMover for spanfile
		AddMembraneMoverOP addmem( new AddMembraneMover( "protocols/membrane/3EFF_TR.span" ) );
		addmem->apply(pose);

		// call test function
		core::Vector pose_tm_com_shouldbe( -0.429458, 0.476042, -1.84529 );
		core::Vector pose_tm_com_is( pose_tm_com( pose ) );

		// test
		position_equal_within_delta( pose_tm_com_is, pose_tm_com_shouldbe , 0.001 );

#endif

	} // pose TM COM

	// residue closest to pose TM COM
	void test_rsd_closest_pose_tm_COM() {

		TR << "=========Test residue closest to pose TM center-of-mass" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		using namespace core::conformation::membrane;
		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;

		// read in pose and create topology object
		Pose pose;
		core::import_pose::pose_from_file( pose, "protocols/membrane/3EFF_TR.pdb" , core::import_pose::PDB_file);

		// call AddMembraneMover for spanfile
		AddMembraneMoverOP addmem( new AddMembraneMover( "protocols/membrane/3EFF_TR.span" ) );
		addmem->apply(pose);

		// test
		TS_ASSERT_EQUALS( rsd_closest_to_pose_tm_com( pose ), 54 );

#endif

	} // rsd closest to pose TM COM

	// create membrane foldtree anchor TM COM with root anchor at pose TM COM
	void test_create_membrane_foldtree_anchor_pose_tmcom() {

		TR << "=========Test create membrane foldtree anchor at TM center-of-mass and pose TM COM" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		using namespace core::conformation::membrane;
		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;

		// read in pose and create topology object
		Pose pose;
		core::import_pose::pose_from_file( pose, "protocols/membrane/3EFF_TR.pdb" , core::import_pose::PDB_file);

		// call AddMembraneMover for spanfile
		AddMembraneMoverOP addmem( new AddMembraneMover( "protocols/membrane/3EFF_TR.span" ) );
		addmem->apply(pose);

		// create foldtree
		create_membrane_foldtree_anchor_pose_tmcom( pose );

		// test
		TS_ASSERT_EQUALS( pose.fold_tree().upstream_jump_residue( 1 ), 557 );
		TS_ASSERT_EQUALS( pose.fold_tree().downstream_jump_residue( 1 ), 54 );
		TS_ASSERT_EQUALS( pose.fold_tree().upstream_jump_residue( 2 ), 54 );
		TS_ASSERT_EQUALS( pose.fold_tree().downstream_jump_residue( 2 ), 158 );
		TS_ASSERT_EQUALS( pose.fold_tree().upstream_jump_residue( 3 ), 54 );
		TS_ASSERT_EQUALS( pose.fold_tree().downstream_jump_residue( 3 ), 352 );
		TS_ASSERT_EQUALS( pose.fold_tree().upstream_jump_residue( 4 ), 54 );
		TS_ASSERT_EQUALS( pose.fold_tree().downstream_jump_residue( 4 ), 491 );

#endif

	} // create membrane foldtree anchor TM COM


	// create membrane foldtree anchor TM COM with root anchor at pose TM COM
	void test_create_membrane_docking_foldtree_from_partners() {

		TR << "=========Test create membrane docking foldtree from partners" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		using namespace core::conformation::membrane;
		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;

		// read in pose and create topology object
		Pose pose;
		core::import_pose::pose_from_file( pose, "protocols/membrane/3EFF_TR.pdb" , core::import_pose::PDB_file);

		// call AddMembraneMover for spanfile
		AddMembraneMoverOP addmem( new AddMembraneMover( "protocols/membrane/3EFF_TR.span" ) );
		addmem->apply(pose);

		// create foldtree
		core::Size i_jump = create_membrane_docking_foldtree_from_partners( pose, core::pose::DockingPartners::docking_partners_from_string("KL_MN") );

		// test
		TS_ASSERT_EQUALS( pose.fold_tree().upstream_jump_residue( 1 ), 557 );
		TS_ASSERT_EQUALS( pose.fold_tree().downstream_jump_residue( 1 ), 19 );
		TS_ASSERT_EQUALS( pose.fold_tree().upstream_jump_residue( 2 ), 19 );
		TS_ASSERT_EQUALS( pose.fold_tree().downstream_jump_residue( 2 ), 158 );
		TS_ASSERT_EQUALS( pose.fold_tree().upstream_jump_residue( 3 ), 352 );
		TS_ASSERT_EQUALS( pose.fold_tree().downstream_jump_residue( 3 ), 491 );
		TS_ASSERT_EQUALS( pose.fold_tree().upstream_jump_residue( 4 ), 19 );
		TS_ASSERT_EQUALS( pose.fold_tree().downstream_jump_residue( 4 ), 352 );
		TS_ASSERT_EQUALS( i_jump, 4 );

#endif

	} // create membrane docking foldtree from partners

	// create membrane foldtree anchor TM COM with root anchor at pose TM COM
	void test_create_specific_membrane_foldtree() {

		TR << "=========Test create specific membrane foldtree from anchors" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		using namespace core::conformation::membrane;
		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;

		// read in pose and create topology object
		Pose pose;
		core::import_pose::pose_from_file( pose, "protocols/membrane/1H2S__tr.pdb" , core::import_pose::PDB_file);

		// call AddMembraneMover for spanfile
		AddMembraneMoverOP addmem( new AddMembraneMover( "protocols/membrane/1H2S__tr.span" ) );
		addmem->apply(pose);

		// create anchors vector
		utility::vector1< core::Vector > anchors;
		core::Vector anchor1( 76, 569, 225 );
		core::Vector anchor2( 76, 235, 284 );
		core::Vector anchor3( 235, 360, 509 );
		core::Vector anchor4( 76, 519, 568 );
		anchors.push_back( anchor1 );
		anchors.push_back( anchor2 );
		anchors.push_back( anchor3 );
		anchors.push_back( anchor4 );

		// create foldtree
		core::Size anchor = create_specific_membrane_foldtree( pose, anchors );

		// test
		TS_ASSERT_EQUALS( pose.fold_tree().upstream_jump_residue( 1 ), 569 );
		TS_ASSERT_EQUALS( pose.fold_tree().downstream_jump_residue( 1 ), 76 );
		TS_ASSERT_EQUALS( pose.fold_tree().upstream_jump_residue( 2 ), 76 );
		TS_ASSERT_EQUALS( pose.fold_tree().downstream_jump_residue( 2 ), 235 );
		TS_ASSERT_EQUALS( pose.fold_tree().upstream_jump_residue( 3 ), 235 );
		TS_ASSERT_EQUALS( pose.fold_tree().downstream_jump_residue( 3 ), 360 );
		TS_ASSERT_EQUALS( pose.fold_tree().upstream_jump_residue( 4 ), 76 );
		TS_ASSERT_EQUALS( pose.fold_tree().downstream_jump_residue( 4 ), 519 );
		TS_ASSERT_EQUALS( pose.conformation().membrane_info()->membrane_jump(), 1 );
		TS_ASSERT_EQUALS( anchor, 76 );

#endif
	} // create specific membrane foldtree

	// create_membrane_multi_partner_foldtree_anchor_tmcom
	void test_create_membrane_multi_partner_foldtree_anchor_tmcom() {

		TR << "=========Test create membrane multi-partner foldtree" << std::endl;

		TS_ASSERT( true );

#ifndef  SERIALIZATION

		using namespace core::conformation::membrane;
		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;

		// read in pose and create topology object
		Pose pose;
		core::import_pose::pose_from_file( pose, "protocols/membrane/1H2S__tr.pdb" , core::import_pose::PDB_file);

		// call AddMembraneMover for spanfile
		AddMembraneMoverOP addmem( new AddMembraneMover( "protocols/membrane/1H2S__tr.span" ) );
		addmem->apply(pose);

		// create foldtree
		utility::vector1< utility::vector1< std::string > > partners{ {"A"}, {"B","C"}, {"D"} };
		utility::vector1< core::Size > interface_jumps = create_membrane_multi_partner_foldtree_anchor_tmcom( pose, partners );

		// test
		TS_ASSERT_EQUALS( pose.fold_tree().upstream_jump_residue( 1 ), 569 );
		TS_ASSERT_EQUALS( pose.fold_tree().downstream_jump_residue( 1 ), 76 );
		TS_ASSERT_EQUALS( pose.fold_tree().upstream_jump_residue( 2 ), 76 );
		TS_ASSERT_EQUALS( pose.fold_tree().downstream_jump_residue( 2 ), 235 );
		TS_ASSERT_EQUALS( pose.fold_tree().upstream_jump_residue( 3 ), 235 );
		TS_ASSERT_EQUALS( pose.fold_tree().downstream_jump_residue( 3 ), 360 );
		TS_ASSERT_EQUALS( pose.fold_tree().upstream_jump_residue( 4 ), 76 );
		TS_ASSERT_EQUALS( pose.fold_tree().downstream_jump_residue( 4 ), 519 );
		TS_ASSERT_EQUALS( pose.conformation().membrane_info()->membrane_jump(), 1 );
		TS_ASSERT_EQUALS( interface_jumps[1], 2 );
		TS_ASSERT_EQUALS( interface_jumps[2], 3 );

#endif

	} // create membrane multi-partner foldtree

private:

	/// @brief Position equal within delta (helper method)
	bool position_equal_within_delta( Vector a, Vector b, Real delta ) {

		TS_ASSERT_DELTA( a.x(), b.x(), delta );
		TS_ASSERT_DELTA( a.y(), b.y(), delta );
		TS_ASSERT_DELTA( a.z(), b.z(), delta );

		return true;
	}

}; // Protocols-level membrane utilities class
