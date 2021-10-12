// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_moves/SimpleThreadingMoverTests
/// @brief  test for SimpleThreadingMover
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ResidueType.hh>

#include <protocols/simple_moves/SimpleThreadingMover.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/helical_bundle/MakeBundle.hh> // For convenience, for making a pose
#include <protocols/helical_bundle/BundleParametrizationCalculator.hh> // For convenience, for making a pose

// Utility Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.simple_moves.SimpleThreadingMoverTests.cxxtest.hh");

// --------------- Test Class --------------- //

class SimpleThreadingMoverTests : public CxxTest::TestSuite {

private:
	core::scoring::ScoreFunctionOP scorefxn_;
	core::pose::Pose ab_pose_aho; //Full PDB
	protocols::antibody::AntibodyInfoOP ab_info_north_aho;

public:

	void setUp() {
		using namespace protocols::antibody;

		core_init();
		core::import_pose::pose_from_file(ab_pose_aho, "protocols/antibody/1bln_AB_aho.pdb");
		ab_info_north_aho = AntibodyInfoOP(new AntibodyInfo(ab_pose_aho, AHO_Scheme, North));

	}

	void tearDown() {
	}

	void test_basic_class() {
		using namespace protocols::antibody;
		using namespace protocols::simple_moves;

		core::Size start = ab_info_north_aho->get_CDR_start(l1, ab_pose_aho);
		core::Size end = ab_info_north_aho->get_CDR_end(l1, ab_pose_aho);

		//L1 length is 16, RSSQSIVHSTGNTYLE.
		std::string seq = "ATATATATATATATAT";
		SimpleThreadingMover threader = SimpleThreadingMover(seq, start);
		TS_ASSERT_THROWS_NOTHING( threader.apply(ab_pose_aho) );
		TS_ASSERT_EQUALS(ab_pose_aho.sequence(start, end), seq)


			//With Options

			std::string seq2 = "ASASASASASASASAS";
		threader.set_pack_rounds(3);
		threader.set_pack_neighbors(true);
		threader.set_neighbor_distance(8);
		threader.set_sequence(seq2, start);
		TS_ASSERT_THROWS_NOTHING( threader.apply( ab_pose_aho ));
		TS_ASSERT_EQUALS(ab_pose_aho.sequence(start, end), seq2);
	}

	void test_oneletter() {
		TR << "Starting SimpleThreadingMoverTests::test_oneletter()" << std::endl;

		using namespace protocols::helical_bundle;
		using namespace protocols::simple_moves;

		MakeBundle makebundle;
		makebundle.set_default_helix_length(30);
		makebundle.set_default_crick_params_file( "alpha_helix_100" );
		makebundle.add_helix();

		{
			core::pose::Pose pose;
			makebundle.apply(pose);
			for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
				TS_ASSERT_EQUALS( pose.residue_type(i).base_name(), "ALA" );
			}

			SimpleThreadingMover threader;
			threader.set_sequence_mode("oneletter");
			threader.set_sequence( "RSTX[DSER]LNE", 3 );
			threader.set_pack_rounds(0);
			TS_ASSERT_THROWS_NOTHING( threader.apply(pose) );
			for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
				if ( i < 3 || i > 9 ) TS_ASSERT_EQUALS( pose.residue_type(i).base_name(), "ALA" );
			}
			TS_ASSERT_EQUALS( pose.residue_type(3).base_name(), "ARG" );
			TS_ASSERT_EQUALS( pose.residue_type(4).base_name(), "SER" );
			TS_ASSERT_EQUALS( pose.residue_type(5).base_name(), "THR" );
			TS_ASSERT_EQUALS( pose.residue_type(6).base_name(), "DSER" );
			TS_ASSERT_EQUALS( pose.residue_type(7).base_name(), "LEU" );
			TS_ASSERT_EQUALS( pose.residue_type(8).base_name(), "ASN" );
			TS_ASSERT_EQUALS( pose.residue_type(9).base_name(), "GLU" );
		}

		{
			core::pose::Pose pose;
			makebundle.apply(pose);
			for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
				TS_ASSERT_EQUALS( pose.residue_type(i).base_name(), "ALA" );
			}

			SimpleThreadingMover threader;
			threader.set_sequence_mode("oneletter");
			threader.set_skip_unknown_mutant(true);
			threader.set_pack_rounds(0);
			threader.set_sequence( "RwTX[001]LNE", 3 );
			TS_ASSERT_THROWS_NOTHING( threader.apply(pose) );
			for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
				if ( i < 3 || i > 9 ) TS_ASSERT_EQUALS( pose.residue_type(i).base_name(), "ALA" );
			}
			TS_ASSERT_EQUALS( pose.residue_type(3).base_name(), "ARG" );
			TS_ASSERT_EQUALS( pose.residue_type(4).base_name(), "ALA" );
			TS_ASSERT_EQUALS( pose.residue_type(5).base_name(), "THR" );
			TS_ASSERT_EQUALS( pose.residue_type(6).base_name(), "001" );
			TS_ASSERT_EQUALS( pose.residue_type(7).base_name(), "LEU" );
			TS_ASSERT_EQUALS( pose.residue_type(8).base_name(), "ASN" );
			TS_ASSERT_EQUALS( pose.residue_type(9).base_name(), "GLU" );
		}

		TR << "Finished SimpleThreadingMoverTests::test_oneletter()" << std::endl;
	}

	void test_threeletter() {
		TR << "Starting SimpleThreadingMoverTests::test_threeletter()" << std::endl;

		using namespace protocols::helical_bundle;
		using namespace protocols::simple_moves;

		MakeBundle makebundle;
		makebundle.set_default_helix_length(30);
		makebundle.set_default_crick_params_file( "alpha_helix_100" );
		makebundle.add_helix();

		{
			core::pose::Pose pose;
			makebundle.apply(pose);
			for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
				TS_ASSERT_EQUALS( pose.residue_type(i).base_name(), "ALA" );
			}

			SimpleThreadingMover threader;
			threader.set_sequence_mode("threeletter");
			threader.set_pack_rounds(0);
			threader.set_sequence( "ARG,SER,THR,DSE,LEU,ASN,GLU", 3 );
			TS_ASSERT_THROWS_NOTHING( threader.apply(pose) );
			for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
				if ( i < 3 || i > 9 ) TS_ASSERT_EQUALS( pose.residue_type(i).base_name(), "ALA" );
			}
			TS_ASSERT_EQUALS( pose.residue_type(3).base_name(), "ARG" );
			TS_ASSERT_EQUALS( pose.residue_type(4).base_name(), "SER" );
			TS_ASSERT_EQUALS( pose.residue_type(5).base_name(), "THR" );
			TS_ASSERT_EQUALS( pose.residue_type(6).base_name(), "DSER" );
			TS_ASSERT_EQUALS( pose.residue_type(7).base_name(), "LEU" );
			TS_ASSERT_EQUALS( pose.residue_type(8).base_name(), "ASN" );
			TS_ASSERT_EQUALS( pose.residue_type(9).base_name(), "GLU" );
		}

		{
			core::pose::Pose pose;
			makebundle.apply(pose);
			for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
				TS_ASSERT_EQUALS( pose.residue_type(i).base_name(), "ALA" );
			}

			SimpleThreadingMover threader;
			threader.set_sequence_mode("threeletter");
			threader.set_pack_rounds(0);
			threader.set_skip_unknown_mutant(true);
			threader.set_sequence( "ARG,fakefakefakefake,THR,001,DSER,ASN,GLU", 3 );
			TS_ASSERT_THROWS_NOTHING( threader.apply(pose) );
			for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
				if ( i < 3 || i > 9 ) TS_ASSERT_EQUALS( pose.residue_type(i).base_name(), "ALA" );
			}
			TS_ASSERT_EQUALS( pose.residue_type(3).base_name(), "ARG" );
			TS_ASSERT_EQUALS( pose.residue_type(4).base_name(), "ALA" );
			TS_ASSERT_EQUALS( pose.residue_type(5).base_name(), "THR" );
			TS_ASSERT_EQUALS( pose.residue_type(6).base_name(), "001" );
			TS_ASSERT_EQUALS( pose.residue_type(7).base_name(), "ALA" ); //DSER isn't DSE
			TS_ASSERT_EQUALS( pose.residue_type(8).base_name(), "ASN" );
			TS_ASSERT_EQUALS( pose.residue_type(9).base_name(), "GLU" );
		}

		TR << "Finished SimpleThreadingMoverTests::test_threeletter()" << std::endl;
	}

	void test_basename() {
		TR << "Starting SimpleThreadingMoverTests::test_basename()" << std::endl;

		using namespace protocols::helical_bundle;
		using namespace protocols::simple_moves;

		MakeBundle makebundle;
		makebundle.set_default_helix_length(30);
		makebundle.set_default_crick_params_file( "alpha_helix_100" );
		makebundle.add_helix();

		{
			core::pose::Pose pose;
			makebundle.apply(pose);
			for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
				TS_ASSERT_EQUALS( pose.residue_type(i).base_name(), "ALA" );
			}

			SimpleThreadingMover threader;
			threader.set_sequence_mode("basename");
			threader.set_pack_rounds(0);
			threader.set_sequence( "ARG,SER,THR,DSER,LEU,ASN,GLU", 3 );
			TS_ASSERT_THROWS_NOTHING( threader.apply(pose) );
			for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
				if ( i < 3 || i > 9 ) TS_ASSERT_EQUALS( pose.residue_type(i).base_name(), "ALA" );
			}
			TS_ASSERT_EQUALS( pose.residue_type(3).base_name(), "ARG" );
			TS_ASSERT_EQUALS( pose.residue_type(4).base_name(), "SER" );
			TS_ASSERT_EQUALS( pose.residue_type(5).base_name(), "THR" );
			TS_ASSERT_EQUALS( pose.residue_type(6).base_name(), "DSER" );
			TS_ASSERT_EQUALS( pose.residue_type(7).base_name(), "LEU" );
			TS_ASSERT_EQUALS( pose.residue_type(8).base_name(), "ASN" );
			TS_ASSERT_EQUALS( pose.residue_type(9).base_name(), "GLU" );
		}

		{
			core::pose::Pose pose;
			makebundle.apply(pose);
			for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
				TS_ASSERT_EQUALS( pose.residue_type(i).base_name(), "ALA" );
			}

			SimpleThreadingMover threader;
			threader.set_sequence_mode("basename");
			threader.set_pack_rounds(0);
			threader.set_skip_unknown_mutant(true);
			threader.set_sequence( "ARG,fakefakefakefake,THR,001,DSER,ALSOFAKE,GLU", 3 );
			TS_ASSERT_THROWS_NOTHING( threader.apply(pose) );
			TS_ASSERT_EQUALS( pose.residue_type(3).base_name(), "ARG" );
			TS_ASSERT_EQUALS( pose.residue_type(4).base_name(), "ALA" );
			TS_ASSERT_EQUALS( pose.residue_type(5).base_name(), "THR" );
			TS_ASSERT_EQUALS( pose.residue_type(6).base_name(), "001" );
			TS_ASSERT_EQUALS( pose.residue_type(7).base_name(), "DSER" );
			TS_ASSERT_EQUALS( pose.residue_type(8).base_name(), "ALA" );
			TS_ASSERT_EQUALS( pose.residue_type(9).base_name(), "GLU" );
		}

		TR << "Finished SimpleThreadingMoverTests::test_basename()" << std::endl;
	}

	void test_fullname() {
		TR << "Starting SimpleThreadingMoverTests::test_fullname()" << std::endl;

		using namespace protocols::helical_bundle;
		using namespace protocols::simple_moves;

		MakeBundle makebundle;
		makebundle.set_default_helix_length(30);
		makebundle.set_default_crick_params_file( "alpha_helix_100" );
		makebundle.add_helix();

		{
			core::pose::Pose pose;
			makebundle.apply(pose);
			for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
				TS_ASSERT_EQUALS( pose.residue_type(i).base_name(), "ALA" );
			}

			SimpleThreadingMover threader;
			threader.set_sequence_mode("fullname");
			threader.set_pack_rounds(0);
			threader.set_sequence( "ARG,SER:N_Methylation,THR,DSER,LEU,ASN,GLU", 3 );
			TS_ASSERT_THROWS_NOTHING( threader.apply(pose) );
			for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
				if ( i < 3 || i > 9 ) TS_ASSERT_EQUALS( pose.residue_type(i).base_name(), "ALA" );
			}
			TS_ASSERT_EQUALS( pose.residue_type(3).base_name(), "ARG" );
			TS_ASSERT_EQUALS( pose.residue_type(4).base_name(), "SER" );
			TS_ASSERT_EQUALS( pose.residue_type(4).name(), "SER:N_Methylation" );
			TS_ASSERT_EQUALS( pose.residue_type(5).base_name(), "THR" );
			TS_ASSERT_EQUALS( pose.residue_type(6).base_name(), "DSER" );
			TS_ASSERT_EQUALS( pose.residue_type(7).base_name(), "LEU" );
			TS_ASSERT_EQUALS( pose.residue_type(8).base_name(), "ASN" );
			TS_ASSERT_EQUALS( pose.residue_type(9).base_name(), "GLU" );

			for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
				if ( i!=4 ) {
					TS_ASSERT( pose.residue_type(i).is_n_methylated() == false );
				} else {
					TS_ASSERT( pose.residue_type(i).is_n_methylated() == true );
				}
			}
		}

		{
			core::pose::Pose pose;
			makebundle.apply(pose);
			for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
				TS_ASSERT_EQUALS( pose.residue_type(i).base_name(), "ALA" );
			}

			SimpleThreadingMover threader;
			threader.set_sequence_mode("fullname");
			threader.set_pack_rounds(0);
			threader.set_skip_unknown_mutant(true);
			threader.set_sequence( "ARG,fakefakefakefake,THR:N_Methylation,001,DSER,ALSOFAKE,GLU", 3 );
			TS_ASSERT_THROWS_NOTHING( threader.apply(pose) );
			for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
				if ( i < 3 || i > 9 ) TS_ASSERT_EQUALS( pose.residue_type(i).base_name(), "ALA" );
			}
			TS_ASSERT_EQUALS( pose.residue_type(3).base_name(), "ARG" );
			TS_ASSERT_EQUALS( pose.residue_type(4).base_name(), "ALA" );
			TS_ASSERT_EQUALS( pose.residue_type(5).base_name(), "THR" );
			TS_ASSERT_EQUALS( pose.residue_type(5).name(), "THR:N_Methylation" );
			TS_ASSERT_EQUALS( pose.residue_type(6).base_name(), "001" );
			TS_ASSERT_EQUALS( pose.residue_type(7).base_name(), "DSER" );
			TS_ASSERT_EQUALS( pose.residue_type(8).base_name(), "ALA" );
			TS_ASSERT_EQUALS( pose.residue_type(9).base_name(), "GLU" );

			for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
				if ( i!=5 ) {
					TS_ASSERT( pose.residue_type(i).is_n_methylated() == false );
				} else {
					TS_ASSERT( pose.residue_type(i).is_n_methylated() == true );
				}
			}
		}

		TR << "Finished SimpleThreadingMoverTests::test_fullname()" << std::endl;
	}

};

