// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/membrane/MPLipidAccessibility.cxxtest.hh
/// @brief  High-level unit test for concave shell calculation
/// @author  JKLeman (julia.koehler1982@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Package Headers
#include <protocols/membrane/MPLipidAccessibility.hh>
#include <protocols/membrane/AddMembraneMover.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/import_pose/import_pose.hh>
#include <core/types.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

static basic::Tracer TR("protocols.membrane.MPLipidAccessibility.cxxtest");

/// @brief Unit Test suite for protocols-level membrane util class
class MPLipidAccessibilityTest : public CxxTest::TestSuite {

public: // test functions

	/// @brief Setup the unit test
	void setUp() {
		core_init();
	}

	/// @brief Tear down the unit test
	void tearDown() {}

	// test MP lipid accessibility
	void test_mp_lipid_accessibility() {

		using namespace core;
		using namespace core::conformation::membrane;
		using namespace core::pose;
		using namespace core::import_pose;
		using namespace protocols::membrane;

		TR << "test MPLipidAccessibility" << std::endl;

		// GPCR
		Pose pose1 = Pose();
		pose_from_file( pose1, "protocols/membrane/1U19__tr.pdb" , core::import_pose::PDB_file );
		std::string spanfile1 = "protocols/membrane/1U19__tr.span";

		AddMembraneMoverOP mem1( new AddMembraneMover( spanfile1 ) );
		MPLipidAccessibilityOP lipid1( new MPLipidAccessibility() );
		mem1->apply( pose1 );
		lipid1->apply( pose1 );

		TS_ASSERT_EQUALS( pose1.pdb_info()->bfactor( 157, 2 ), 50.0 );
		TS_ASSERT_EQUALS( pose1.pdb_info()->bfactor( 260, 2 ), 50.0 );
		TS_ASSERT_EQUALS( pose1.pdb_info()->bfactor(  72, 2 ), 0.0 );
		TS_ASSERT_EQUALS( pose1.pdb_info()->bfactor( 127, 2 ), 0.0 );
		TS_ASSERT_EQUALS( pose1.pdb_info()->bfactor( 125, 2 ), 0.0 );

		// ion channel
		Pose pose2 = Pose();
		pose_from_file( pose2, "protocols/membrane/4DXW__tr.pdb" , core::import_pose::PDB_file );
		std::string spanfile2 = "protocols/membrane/4DXW__tr.span";

		AddMembraneMoverOP mem2( new AddMembraneMover( spanfile2 ) );
		MPLipidAccessibilityOP lipid2( new MPLipidAccessibility() );
		mem2->apply( pose2 );
		lipid2->apply( pose2 );

		TS_ASSERT_EQUALS( pose2.pdb_info()->bfactor( 108, 2 ), 50.0 );
		TS_ASSERT_EQUALS( pose2.pdb_info()->bfactor( 780, 2 ), 50.0 );
		TS_ASSERT_EQUALS( pose2.pdb_info()->bfactor( 575, 2 ), 0.0 );
		TS_ASSERT_EQUALS( pose2.pdb_info()->bfactor( 413, 2 ), 0.0 );
		TS_ASSERT_EQUALS( pose2.pdb_info()->bfactor( 311, 2 ), 0.0 );

		// chloride channel
		Pose pose3 = Pose();
		pose_from_file( pose3, "protocols/membrane/4ENE__tr.pdb" , core::import_pose::PDB_file );
		std::string spanfile3 = "protocols/membrane/4ENE__tr.span";

		AddMembraneMoverOP mem3( new AddMembraneMover( spanfile3 ) );
		MPLipidAccessibilityOP lipid3( new MPLipidAccessibility() );
		mem3->apply( pose3 );
		lipid3->apply( pose3 );

		TS_ASSERT_EQUALS( pose3.pdb_info()->bfactor( 203, 2 ), 50.0 );
		TS_ASSERT_EQUALS( pose3.pdb_info()->bfactor( 141, 2 ), 50.0 );
		TS_ASSERT_EQUALS( pose3.pdb_info()->bfactor( 120, 2 ), 50.0 );
		TS_ASSERT_EQUALS( pose3.pdb_info()->bfactor( 582, 2 ), 0.0 );
		TS_ASSERT_EQUALS( pose3.pdb_info()->bfactor( 347, 2 ), 0.0 );

		// ABC transporter
		Pose pose4 = Pose();
		pose_from_file( pose4, "protocols/membrane/3TUI__tr.pdb" , core::import_pose::PDB_file );
		std::string spanfile4 = "protocols/membrane/3TUI__tr.span";

		AddMembraneMoverOP mem4( new AddMembraneMover( spanfile4 ) );
		MPLipidAccessibilityOP lipid4( new MPLipidAccessibility() );
		mem4->apply( pose4 );
		lipid4->apply( pose4 );

		TS_ASSERT_EQUALS( pose4.pdb_info()->bfactor( 148, 2 ), 50.0 );
		TS_ASSERT_EQUALS( pose4.pdb_info()->bfactor( 319, 2 ), 50.0 );
		TS_ASSERT_EQUALS( pose4.pdb_info()->bfactor( 368, 2 ), 0.0 );
		TS_ASSERT_EQUALS( pose4.pdb_info()->bfactor( 103, 2 ), 0.0 );
		TS_ASSERT_EQUALS( pose4.pdb_info()->bfactor( 374, 2 ), 0.0 );

		// beta-barrel
		Pose pose5 = Pose();
		pose_from_file( pose5, "protocols/membrane/3CSL__tr.pdb" , core::import_pose::PDB_file );
		std::string spanfile5 = "protocols/membrane/3CSL__tr.span";

		AddMembraneMoverOP mem5( new AddMembraneMover( spanfile5 ) );
		MPLipidAccessibilityOP lipid5( new MPLipidAccessibility() );
		mem5->apply( pose5 );
		lipid5->apply( pose5 );

		TS_ASSERT_EQUALS( pose5.pdb_info()->bfactor( 753, 2 ), 50.0 );
		TS_ASSERT_EQUALS( pose5.pdb_info()->bfactor( 667, 2 ), 50.0 );
		TS_ASSERT_EQUALS( pose5.pdb_info()->bfactor( 720, 2 ), 50.0 );
		TS_ASSERT_EQUALS( pose5.pdb_info()->bfactor( 392, 2 ), 0.0 );
		TS_ASSERT_EQUALS( pose5.pdb_info()->bfactor( 5, 2 ), 0.0 );

	}

	/// @brief test parse_my_tag() read-in
	void test_parse_my_tag_input() {
		using namespace protocols::membrane;

		std::stringstream tag_ss("<MPLipidAccessibility angle_cutoff=64.0 slice_width=9.0 shell_radius=5.0 dist_cutoff=9.0 tm_alpha=false />");
		utility::tag::TagCOP tag = utility::tag::Tag::create( tag_ss );

		MPLipidAccessibilityOP xlip( new MPLipidAccessibility() );
		xlip->parse_my_tag( tag );

		TS_ASSERT_EQUALS(xlip->get_angle_cutoff(), 64.0);
		TS_ASSERT_EQUALS(xlip->get_slice_width(), 9.0);
		TS_ASSERT_EQUALS(xlip->get_shell_radius(), 5.0);
		TS_ASSERT_EQUALS(xlip->get_dist_cutoff(), 9.0);
		TS_ASSERT_EQUALS(xlip->get_tm_alpha(), false);
	}

};
