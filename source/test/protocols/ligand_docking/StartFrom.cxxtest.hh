// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/pose/util.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <protocols/ligand_docking/StartFrom.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>


static basic::Tracer TR("protocols.ligand_docking.StartFrom.cxxtest");


class StartFromTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
		// Residue definitions can't be supplied on the command line b/c
		// the ResidueTypeSet is already initialized.
		using namespace core::chemical;
		utility::vector1< std::string > params_files;
		ResidueTypeSetCOP const_residue_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		ResidueTypeSet & residue_set = const_cast< ResidueTypeSet & >(*const_residue_set);
		if(!residue_set.has_name("ZNx")) params_files.push_back("protocols/ligand_docking/ZNx.params");
		if(!residue_set.has_name("CP1")) params_files.push_back("protocols/ligand_docking/7cpa.params");
		residue_set.read_files(params_files);
	}

	void tearDown() {}

	void test_specified_position() {
		using namespace protocols::ligand_docking;

		core::pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, "protocols/ligand_docking/7cpa_7cpa_native.pdb" );

		StartFrom mover;
		mover.chain("X");

		core::Vector v(1.0,20.0,100.0);
		mover.add_coords(v); //default

		mover.apply(pose);

		core::Size jump_id = core::pose::get_jump_id_from_chain("X", pose);
		core::Vector out( protocols::geometry::downstream_centroid_by_jump(pose, jump_id) );
		TS_ASSERT_DELTA( out.x(), v.x(), 1e-4 );
		TS_ASSERT_DELTA( out.y(), v.y(), 1e-4 );
		TS_ASSERT_DELTA( out.z(), v.z(), 1e-4 );
	}

	void test_specified_position_nbr() {
		using namespace protocols::ligand_docking;

		core::pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, "protocols/ligand_docking/7cpa_7cpa_native.pdb" );

		StartFrom mover;
		mover.chain("X");
		mover.use_nbr( true );

		core::Vector v(1.0,20.0,100.0);
		mover.add_coords(v); //default

		mover.apply(pose);

		core::Vector out( pose.residue(309).xyz("N1") );
		TS_ASSERT_DELTA( out.x(), v.x(), 1e-4 );
		TS_ASSERT_DELTA( out.y(), v.y(), 1e-4 );
		TS_ASSERT_DELTA( out.z(), v.z(), 1e-4 );
	}

	void test_multiple_positions() {
		using namespace protocols::ligand_docking;

		core::pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, "protocols/ligand_docking/7cpa_7cpa_native.pdb" );

		StartFrom mover;
		mover.chain("X");
		mover.use_nbr( true );

		for( core::Real pp(0); pp <= 25; pp += 5 ) {
			core::Vector v(pp,20.0,100.0);
			mover.add_coords(v); //default
		}

		std::map< core::Size, core::Size > counts;
		for( core::Size ii(1); ii <= 200; ++ii ) {
			mover.apply(pose);
			core::Vector out( pose.residue(309).xyz("N1") );
			// Deliberate conversion/rounding to eliminate noise
			core::Size x_pos( pose.residue(309).xyz("N1").x()+0.25 );
			++counts[ x_pos ];
		}

		TS_ASSERT_EQUALS( counts.size(), 6 ); // Only the specified coordinates
		//Each should get roughly 200/6 = 33 items
		for( core::Size pp(0); pp <= 25; pp += 5 ) {
			TS_ASSERT_LESS_THAN( 23, counts[pp] );
			TS_ASSERT_LESS_THAN( counts[pp], 43 );
		}
	}

	// Test that non-applicable pdb tags or hashes don't change the items
	void test_spurious_items() {
		using namespace protocols::ligand_docking;

		core::pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, "protocols/ligand_docking/7cpa_7cpa_native.pdb" );

		StartFrom mover;
		mover.chain("X");

		core::Vector v(1.0,20.0,100.0);
		mover.add_coords(v); //default

		// Some items which shouldn't be used
		core::Vector pv( 99.0, 96.0, 33.0 );
		mover.add_coords(pv, "some_other.pdb");
		core::Vector hv( 199.0, -96.0, 333.0 );
		mover.add_coords_hash(hv, "aa2aff055d19bc32e483df7ff4ae08361a768931");

		for( core::Size ii(1); ii <= 10; ++ii ) {
			mover.apply(pose);
			core::Size jump_id = core::pose::get_jump_id_from_chain("X", pose);
			core::Vector out( protocols::geometry::downstream_centroid_by_jump(pose, jump_id) );
			// Should get the single default one
			TS_ASSERT_DELTA( out.x(), v.x(), 1e-4 );
			TS_ASSERT_DELTA( out.y(), v.y(), 1e-4 );
			TS_ASSERT_DELTA( out.z(), v.z(), 1e-4 );
		}
	}

	void test_pdb_tag() {
		using namespace protocols::ligand_docking;

		core::pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, "protocols/ligand_docking/7cpa_7cpa_native.pdb" );

		StartFrom mover;
		mover.chain("X");

		core::Vector v(1.0,20.0,100.0);
		mover.add_coords(v); //default
		core::Vector pv( 99.0, 96.0, 33.0 );
		mover.add_coords(pv, "EMPTY_JOB_use_jd2"); // This is the tag used during unit tests
		core::Vector hv( 199.0, -96.0, 333.0 );
		mover.add_coords_hash(hv, "aa2aff055d19bc32e483df7ff4ae08361a768931");

		for( core::Size ii(1); ii <= 10; ++ii ) {
			mover.apply(pose);
			core::Size jump_id = core::pose::get_jump_id_from_chain("X", pose);
			core::Vector out( protocols::geometry::downstream_centroid_by_jump(pose, jump_id) );
			// Should get the PDB tagged one for preference
			TS_ASSERT_DELTA( out.x(), pv.x(), 1e-4 );
			TS_ASSERT_DELTA( out.y(), pv.y(), 1e-4 );
			TS_ASSERT_DELTA( out.z(), pv.z(), 1e-4 );
		}
	}

	void test_hash_tag() {
		using namespace protocols::ligand_docking;

		core::pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, "protocols/ligand_docking/7cpa_7cpa_native.pdb" );

		StartFrom mover;
		mover.chain("X");

		core::Vector v(1.0,20.0,100.0);
		mover.add_coords(v); //default
		core::Vector pv( 99.0, 96.0, 33.0 );
		mover.add_coords(pv, "some_other.pdb");
		core::Vector hv( 199.0, -96.0, 333.0 );
		mover.add_coords_hash(hv, core::pose::get_sha1_hash_excluding_chain('X',pose) );

		for( core::Size ii(1); ii <= 10; ++ii ) {
			mover.apply(pose);
			core::Size jump_id = core::pose::get_jump_id_from_chain("X", pose);
			core::Vector out( protocols::geometry::downstream_centroid_by_jump(pose, jump_id) );
			// Should get the one from the hash for preference
			TS_ASSERT_DELTA( out.x(), hv.x(), 1e-4 );
			TS_ASSERT_DELTA( out.y(), hv.y(), 1e-4 );
			TS_ASSERT_DELTA( out.z(), hv.z(), 1e-4 );
		}
	}

	// Test the ability to read coordinates from a PDB file
	void test_PDB_file() {
		using namespace protocols::ligand_docking;

		core::pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, "protocols/ligand_docking/7cpa_7cpa_native.pdb" );

		StartFrom mover;
		mover.chain("X");
		mover.parse_pdb_file("protocols/ligand_docking/ligand_pos.pdb");

		std::map< core::Size, core::Size > counts;
		for( core::Size ii(1); ii <= 100; ++ii ) {
			mover.apply(pose);
			core::Size jump_id = core::pose::get_jump_id_from_chain("X", pose);
			core::Vector out( protocols::geometry::downstream_centroid_by_jump(pose, jump_id) );
			// Deliberate conversion/rounding to eliminate noise
			core::Size y_pos( out.y()+0.25 );
			++counts[ y_pos ];
		}

		utility::vector1< core::Size> y_pos;
		y_pos.push_back( 8 );
		y_pos.push_back( 26 );
		y_pos.push_back( 27 );
		y_pos.push_back( 35 );
		//y_pos.push_back( 38 ); // Hydrogen - skipped
		//y_pos.push_back( 39 ); // Hydrogen - skipped

		TS_ASSERT_EQUALS( counts.size(), y_pos.size() ); // Only the specified coordinates
		//Each should get roughly 100/4 = 25 items
		for( core::Size pp(1); pp <= y_pos.size(); ++pp ) {
			TS_ASSERT_LESS_THAN( 17, counts[y_pos[pp]] );
			TS_ASSERT_LESS_THAN( counts[y_pos[pp]], 33 );
		}
	}

	// Test the ability to read coordinates from a PDB file for only a particular atom
	void test_PDB_file_atom_name() {
		using namespace protocols::ligand_docking;

		core::pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, "protocols/ligand_docking/7cpa_7cpa_native.pdb" );

		StartFrom mover;
		mover.chain("X");
		mover.parse_pdb_file("protocols/ligand_docking/ligand_pos.pdb", "O1");

		std::map< core::Size, core::Size > counts;
		for( core::Size ii(1); ii <= 50; ++ii ) {
			mover.apply(pose);
			core::Size jump_id = core::pose::get_jump_id_from_chain("X", pose);
			core::Vector out( protocols::geometry::downstream_centroid_by_jump(pose, jump_id) );
			// Deliberate conversion/rounding to eliminate noise
			core::Size y_pos( out.y()+0.25 );
			++counts[ y_pos ];
		}

		utility::vector1< core::Size> y_pos;
		y_pos.push_back( 8 );
		//y_pos.push_back( 26 );
		//y_pos.push_back( 27 );
		y_pos.push_back( 35 );
		//y_pos.push_back( 38 ); // Hydrogen - skipped
		//y_pos.push_back( 39 ); // Hydrogen - skipped

		TS_ASSERT_EQUALS( counts.size(), y_pos.size() ); // Only the specified coordinates
		//Each should get roughly 50/2 = 25 items
		for( core::Size pp(1); pp <= y_pos.size(); ++pp ) {
			TS_ASSERT_LESS_THAN( 17, counts[y_pos[pp]] );
			TS_ASSERT_LESS_THAN( counts[y_pos[pp]], 33 );
		}
	}

	// Test the ability to read coordinates from a json file
	void test_json_file() {
		using namespace protocols::ligand_docking;

		core::pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, "protocols/ligand_docking/7cpa_7cpa_native.pdb" );

		StartFrom mover;
		mover.chain("X");
		mover.parse_startfrom_file("protocols/ligand_docking/ligand_pos.json");

		std::map< core::Size, core::Size > counts;
		for( core::Size ii(1); ii <= 50; ++ii ) {
			mover.apply(pose);
			core::Size jump_id = core::pose::get_jump_id_from_chain("X", pose);
			core::Vector out( protocols::geometry::downstream_centroid_by_jump(pose, jump_id) );
			// Deliberate conversion/rounding to eliminate noise
			core::Size x_pos( out.x()+0.25 );
			++counts[ x_pos ];
		}

		utility::vector1< core::Size> x_pos;
		x_pos.push_back( 10 ); // file_name works
		x_pos.push_back( 20 ); // input_tag works
		//x_pos.push_back( 30 ); // hash value overrules file_name/input_tag
		//x_pos.push_back( 40 ); // hash only doesn't count
		//x_pos.push_back( 50 ); // nor does the default
		//x_pos.push_back( 60 ); // or another tag

		TS_ASSERT_EQUALS( counts.size(), x_pos.size() ); // Only the specified coordinates
		//Each should get roughly 50/2 = 25 items
		for( core::Size pp(1); pp <= x_pos.size(); ++pp ) {
			TS_ASSERT_LESS_THAN( 17, counts[x_pos[pp]] );
			TS_ASSERT_LESS_THAN( counts[x_pos[pp]], 33 );
		}
	}

};

