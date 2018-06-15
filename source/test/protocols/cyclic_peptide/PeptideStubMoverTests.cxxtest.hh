// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/cyclic_peptide/PeptideStubMoverTests.cxxtest.hh
/// @brief  Unit tests for the PeptideStubMover, used to build peptide geometry.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/util/rosettascripts.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/cyclic_peptide/PeptideStubMover.hh>

// Protocols Headers
#include <protocols/simple_moves/MutateResidue.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <utility/pointer/memory.hh>

static basic::Tracer TR("PeptideStubMoverTests");


class PeptideStubMoverTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options( "-write_all_connect_info true" );
		core::pose::PoseOP pose( utility::pointer::make_shared< core::pose::Pose >() );
		core::import_pose::pose_from_file( *pose, "protocols/cyclic_peptide/peptidestubmover_test_pose.pdb", false, core::import_pose::PDB_file );
		for ( core::Size i(1), imax(pose->total_residue()); i<=imax; ++i ) {
			protocols::simple_moves::MutateResidue mutres( i, "VAL" ); //Mutate the whole pose to valine.
			mutres.apply(*pose);
		}
		core::pose::add_lower_terminus_type_to_pose_residue(*pose, 1);
		core::pose::add_upper_terminus_type_to_pose_residue(*pose, pose->size());

		pose->update_residue_neighbors();
		pose_ = pose; //Nonconst to const.
	}

	void tearDown(){

	}

	void test_append() {
		core::pose::Pose pose( *(pose_) ); //Local copy

		core::Size const original_pose_size( pose.total_residue() );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Append", "ALA", 0, false, "", 0, pose.total_residue(), "" );
		stubmover.apply(pose);
		TS_ASSERT_EQUALS(pose.total_residue(), original_pose_size+1 );
		TS_ASSERT_EQUALS( pose.residue_type( pose.total_residue() ).name3(), "ALA" );
		TS_ASSERT_EQUALS( pose.residue_type( original_pose_size ).name3(), "VAL" );
		TS_ASSERT_EQUALS( pose.residue(pose.total_residue()).connected_residue_at_lower(), original_pose_size );
		TS_ASSERT_EQUALS( pose.residue(original_pose_size).connected_residue_at_upper(), pose.total_residue() );
		TS_ASSERT_DELTA( pose.residue( original_pose_size ).xyz("C").distance( pose.residue( original_pose_size + 1 ).xyz("N") ), 1.328685, 0.01 ); //Check peptide bond length.
	}

	void test_append2() {
		core::pose::Pose pose( *(pose_) ); //Local copy

		core::Size const original_pose_size( pose.total_residue() );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Append", "ALA", 0, false, "", 0, 0, "" ); //Should automatically select last
		stubmover.apply(pose);
		TS_ASSERT_EQUALS(pose.total_residue(), original_pose_size+1 );
		TS_ASSERT_EQUALS( pose.residue_type( pose.total_residue() ).name3(), "ALA" );
		TS_ASSERT_EQUALS( pose.residue_type( original_pose_size ).name3(), "VAL" );
		TS_ASSERT_EQUALS( pose.residue(pose.total_residue()).connected_residue_at_lower(), original_pose_size );
		TS_ASSERT_EQUALS( pose.residue(original_pose_size).connected_residue_at_upper(), pose.total_residue() );
		TS_ASSERT_DELTA( pose.residue( original_pose_size ).xyz("C").distance( pose.residue( original_pose_size + 1 ).xyz("N") ), 1.328685, 0.01 ); //Check peptide bond length.
	}

	void test_multi_append() {
		core::pose::Pose pose( *(pose_) ); //Local copy

		core::Size const original_pose_size( pose.total_residue() );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		for ( core::Size i(0); i<5; ++i ) {
			stubmover.add_residue( "Append", "ALA", 0, false, "", 0, pose.total_residue() + i, "" );
		}
		stubmover.apply(pose);
		TS_ASSERT_EQUALS(pose.total_residue(), original_pose_size+5 );
		TS_ASSERT_EQUALS( pose.residue_type( original_pose_size ).name3(), "VAL" );
		for ( core::Size i(1); i<=5; ++i ) {
			TS_ASSERT_EQUALS( pose.residue_type( original_pose_size + i ).name3(), "ALA" );
			TS_ASSERT_EQUALS( pose.residue( original_pose_size + i ).connected_residue_at_lower(), original_pose_size + i - 1 );
			TS_ASSERT_EQUALS( pose.residue( original_pose_size + i - 1 ).connected_residue_at_upper(), original_pose_size + i );
			TS_ASSERT_DELTA( pose.residue( original_pose_size + i - 1 ).xyz("C").distance( pose.residue( original_pose_size + i ).xyz("N") ), 1.328685, 0.01 ); //Check peptide bond length.
		}
	}

	void test_append_repeat() {
		core::pose::Pose pose( *(pose_) ); //Local copy

		core::Size const original_pose_size( pose.total_residue() );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Append", "ALA", 0, false, "", 5, pose.total_residue(), "" );
		stubmover.apply(pose);
		TS_ASSERT_EQUALS(pose.total_residue(), original_pose_size+5 );
		TS_ASSERT_EQUALS( pose.residue_type( original_pose_size ).name3(), "VAL" );
		for ( core::Size i(1); i<=5; ++i ) {
			TS_ASSERT_EQUALS( pose.residue_type( original_pose_size + i ).name3(), "ALA" );
			TS_ASSERT_EQUALS( pose.residue( original_pose_size + i ).connected_residue_at_lower(), original_pose_size + i - 1 );
			TS_ASSERT_EQUALS( pose.residue( original_pose_size + i - 1 ).connected_residue_at_upper(), original_pose_size + i );
			TS_ASSERT_DELTA( pose.residue( original_pose_size + i - 1 ).xyz("C").distance( pose.residue( original_pose_size + i ).xyz("N") ), 1.328685, 0.01 ); //Check peptide bond length.
		}
	}

	void test_prepend() {
		core::pose::Pose pose( *(pose_) ); //Local copy

		core::Size const original_pose_size( pose.total_residue() );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Prepend", "ALA", 0, false, "", 0, 1, "" );
		stubmover.apply(pose);
		TS_ASSERT_EQUALS(pose.total_residue(), original_pose_size+1 );
		TS_ASSERT_EQUALS( pose.residue_type( pose.total_residue() ).name3(), "VAL" );
		TS_ASSERT_EQUALS( pose.residue_type( original_pose_size ).name3(), "VAL" );
		TS_ASSERT_EQUALS( pose.residue_type( 2 ).name3(), "VAL" );
		TS_ASSERT_EQUALS( pose.residue_type( 1 ).name3(), "ALA" );
		TS_ASSERT_EQUALS( pose.residue(2).connected_residue_at_lower(), 1 );
		TS_ASSERT_EQUALS( pose.residue(1).connected_residue_at_upper(), 2 );
		TS_ASSERT_DELTA( pose.residue(1).xyz("C").distance( pose.residue(2).xyz("N") ), 1.328685, 0.01 ); //Check peptide bond length.
	}

	void test_prepend2() {
		core::pose::Pose pose( *(pose_) ); //Local copy

		core::Size const original_pose_size( pose.total_residue() );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Prepend", "ALA", 0, false, "", 0, 0, "" ); //Should auto-detect?
		stubmover.apply(pose);
		TS_ASSERT_EQUALS(pose.total_residue(), original_pose_size+1 );
		TS_ASSERT_EQUALS( pose.residue_type( pose.total_residue() ).name3(), "VAL" );
		TS_ASSERT_EQUALS( pose.residue_type( original_pose_size ).name3(), "VAL" );
		TS_ASSERT_EQUALS( pose.residue_type( 2 ).name3(), "VAL" );
		TS_ASSERT_EQUALS( pose.residue_type( 1 ).name3(), "ALA" );
		TS_ASSERT_EQUALS( pose.residue(2).connected_residue_at_lower(), 1 );
		TS_ASSERT_EQUALS( pose.residue(1).connected_residue_at_upper(), 2 );
		TS_ASSERT_DELTA( pose.residue(1).xyz("C").distance( pose.residue(2).xyz("N") ), 1.328685, 0.01 ); //Check peptide bond length.
	}

	void test_multi_prepend() {
		core::pose::Pose pose( *(pose_) ); //Local copy

		core::Size const original_pose_size( pose.total_residue() );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		for ( core::Size i(0); i<5; ++i ) {
			stubmover.add_residue( "Prepend", "ALA", 0, false, "", 0, 1, "" );
		}
		stubmover.apply(pose);
		TS_ASSERT_EQUALS(pose.total_residue(), original_pose_size+5 );
		TS_ASSERT_EQUALS( pose.residue_type( pose.total_residue() ).name3(), "VAL" );
		TS_ASSERT_EQUALS( pose.residue_type( 6 ).name3(), "VAL" );
		for ( core::Size i(1); i<=5; ++i ) {
			TS_ASSERT_EQUALS( pose.residue_type( i ).name3(), "ALA" );
			TS_ASSERT_EQUALS( pose.residue( i+1 ).connected_residue_at_lower(), i );
			TS_ASSERT_EQUALS( pose.residue( i ).connected_residue_at_upper(), i+1 );
			TS_ASSERT_DELTA( pose.residue(i).xyz("C").distance( pose.residue(i+1).xyz("N") ), 1.328685, 0.01 ); //Check peptide bond length.
		}
	}

	void test_prepend_repeat() {
		core::pose::Pose pose( *(pose_) ); //Local copy

		core::Size const original_pose_size( pose.total_residue() );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Prepend", "ALA", 0, false, "", 5, 1, "" );
		stubmover.apply(pose);
		TS_ASSERT_EQUALS(pose.total_residue(), original_pose_size+5 );
		TS_ASSERT_EQUALS( pose.residue_type( pose.total_residue() ).name3(), "VAL" );
		TS_ASSERT_EQUALS( pose.residue_type( 6 ).name3(), "VAL" );
		for ( core::Size i(1); i<=5; ++i ) {
			TS_ASSERT_EQUALS( pose.residue_type( i ).name3(), "ALA" );
			TS_ASSERT_EQUALS( pose.residue( i+1 ).connected_residue_at_lower(), i );
			TS_ASSERT_EQUALS( pose.residue( i ).connected_residue_at_upper(), i+1 );

			TS_ASSERT_DELTA( pose.residue(i).xyz("C").distance( pose.residue(i+1).xyz("N") ), 1.328685, 0.01 ); //Check peptide bond length.
		}
	}

	void test_prepend_fail() {
		core::pose::Pose pose( *(pose_) ); //Local copy
		core::Size const original_pose_size( pose.total_residue() );

		basic::datacache::DataMap data;
		Filters_map filters;
		Movers_map movers;
		TagCOP tag = tagptr_from_string("<PeptideStubMover name=\"psm\" reset=\"false\">\n<Append resname=\"ALA\" anchor_rsd=\"72\"/>\n <Prepend resname=\"ALA\" anchor_rsd=\"1\" repeat=\"3\"/></PeptideStubMover>\n");

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.parse_my_tag( tag, data, filters, movers, pose );
		stubmover.apply(pose);
		TS_ASSERT_EQUALS( pose.total_residue(), original_pose_size + 4 );
		TS_ASSERT_EQUALS( pose.residue_type( 1 ).name3(), "ALA" );

	}


private:

	core::pose::PoseCOP pose_;



};
