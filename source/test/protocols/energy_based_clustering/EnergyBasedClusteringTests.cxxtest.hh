// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/energy_based_clustering/EnergyBasedClusteringTests.cxxtest.hh
/// @brief  Unit tests for the EnergyBasedClusteringProtocol.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers


// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/import_pose/import_pose.hh>

// Protocols Headers
#include <protocols/energy_based_clustering/EnergyBasedClusteringProtocol.hh>
#include <protocols/energy_based_clustering/EnergyBasedClusteringOptions.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("EnergyBasedClusteringTests");

class EnergyBasedClusteringTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options( "-energy_based_clustering:use_CB true" );

	}

	void tearDown(){

	}

	/// @brief Test that we're counting asymmetric binstrings correctly.
	void test_num_asymmetric_binstrings() {
		using namespace protocols::energy_based_clustering;
		std::set< std::string > binstrings1;
		binstrings1.insert( "AAXBYYBB" ); //Asymmetric
		binstrings1.insert( "XXXBYYBB" ); //Asymmetric
		binstrings1.insert( "YYBBBBYY" ); //Symmetric
		binstrings1.insert( "AAXBAAXB" ); //Symmetric
		binstrings1.insert( "YYXBBBAX" ); //Asymmetric
		binstrings1.insert( "YYXBBBAY" ); //Asymmetric

		std::set< std::string > binstrings2;
		binstrings2.insert( "ABOZXYZO" ); //Symmetric
		binstrings2.insert( "BBAZYYXO" ); //Symmetric
		binstrings2.insert( "YYBBBBYY" ); //Symmetric
		binstrings2.insert( "XXXBYYBB" ); //Asymmetric
		binstrings2.insert( "AXXXXAAA" ); //Symmetric

		std::set< std::string > binstrings3;
		binstrings3.insert( "ABOZXYZO" ); //Symmetric
		binstrings3.insert( "BBAZYYXO" ); //Symmetric

		TS_ASSERT_EQUALS( binstrings1.size(), 6 );
		TS_ASSERT_EQUALS( binstrings2.size(), 5 );
		TS_ASSERT_EQUALS( binstrings3.size(), 2 );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::num_asymmetric_binstrings( binstrings1, true ), 4 );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::num_asymmetric_binstrings( binstrings2, true ), 1 );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::num_asymmetric_binstrings( binstrings3, true ), 0 );
	}

	/// @brief Test that we're correctly counting the number of asymmetric bin strings that have their mirror
	/// counterpart represented in the set.
	void test_num_asymmetric_binstrings_with_mirror_counterpart_represented() {
		using namespace protocols::energy_based_clustering;
		std::set< std::string > binstrings1;
		binstrings1.insert( "AAXBYYBB" ); //Asymmetric -- paired
		binstrings1.insert( "XAXBYYBB" ); //Asymmetric
		binstrings1.insert( "YAXBYYBB" ); //Asymmetric --paired
		binstrings1.insert( "BBAZYYXO" ); //Symmetric
		binstrings1.insert( "YBBYYBXA" ); //Asymmetric -- paired
		binstrings1.insert( "BYYXXAYB" ); //Asymmetric -- paired
		binstrings1.insert( "AAXBYOBB" ); //Asymmetric
		binstrings1.insert( "AAXBYZBB" ); //Asymmetric
		binstrings1.insert( "ABOZXYZO" ); //Symmetric
		binstrings1.insert( "AXXXXAAA" ); //Symmetric

		std::set< std::string > binstrings2;
		binstrings2.insert( "ABOOXYZZ" ); //Symmetric
		binstrings2.insert( "XAXBYYBB" ); //Asymmetric
		binstrings2.insert( "YAXBYOOB" ); //Asymmetric
		binstrings2.insert( "BBAZYYXO" ); //Symmetric
		binstrings2.insert( "YBBYYBXZ" ); //Asymmetric
		binstrings2.insert( "BYYXXAYB" ); //Asymmetric
		binstrings2.insert( "AAXBYOBB" ); //Asymmetric
		binstrings2.insert( "AAXBYZBB" ); //Asymmetric
		binstrings2.insert( "ABOZXYZO" ); //Symmetric
		binstrings2.insert( "AXXXXAAA" ); //Symmetric

		TS_ASSERT_EQUALS( binstrings1.size(), 10 );
		TS_ASSERT_EQUALS( binstrings2.size(), 10 );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::num_asymmetric_binstrings_with_mirror_counterpart_represented( binstrings1, true ), 4 );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::num_asymmetric_binstrings_with_mirror_counterpart_represented( binstrings2, true ), 0 );
	}

	/// @brief Test that we're assigning ABOXYZ bins correctly.
	void test_determine_ABOXYZ_bin() {
		using namespace protocols::energy_based_clustering;
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::determine_ABOXYZ_bin( -61.1, -41.1, 182.3 ), 'A' );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::determine_ABOXYZ_bin( -61.1, -41.1, 173.5 ), 'A' );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::determine_ABOXYZ_bin( -61.1, -43.1, -173.5 ), 'A' );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::determine_ABOXYZ_bin( 61.1, 41.1, 182.3 ), 'X' );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::determine_ABOXYZ_bin( 61.1, 41.1, 173.5 ), 'X' );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::determine_ABOXYZ_bin( 61.1, 43.1, -173.5 ), 'X' );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::determine_ABOXYZ_bin( -135, 135, 173.5 ), 'B' );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::determine_ABOXYZ_bin( 135-360.0, -135, 173.5 ), 'Y' );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::determine_ABOXYZ_bin( -135+360.0, 50, 173.5 ), 'A' );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::determine_ABOXYZ_bin( -135, -80, 173.5 ), 'B' );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::determine_ABOXYZ_bin( 135, -50, 173.5 ), 'X' );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::determine_ABOXYZ_bin( 135, -80, 173.5 ), 'Y' );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::determine_ABOXYZ_bin( -135, -80, 3.5 ), 'O' );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::determine_ABOXYZ_bin( 135, -80, 3.5 ), 'Z' );

		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::get_mirror_bin('A'), 'X' );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::get_mirror_bin('B'), 'Y' );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::get_mirror_bin('O'), 'Z' );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::get_mirror_bin('X'), 'A' );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::get_mirror_bin('Y'), 'B' );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::get_mirror_bin('Z'), 'O' );
		TR << "The following should throw..." << std::endl;
		TS_ASSERT_THROWS_ANYTHING( EnergyBasedClusteringProtocol::get_mirror_bin('Q') );

		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::get_mirror_bin_sequence( "BBABXYYBAXYBBYOBOZBYA", false ), "YYXYABBYXABYYBZYZOYBX" );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::get_mirror_bin_sequence( "BBABXYYBAXYBBYOBOZBYA", true ), "ABBYXABYYBZYZOYBXYYXY" );
	}

	/// @brief Test that the permute_string function works correctly.
	void test_permute_string() {
		using namespace protocols::energy_based_clustering;

		std::string const string1(""), string2("Q"), string3("Toronto");
		std::string const string1b( EnergyBasedClusteringProtocol::permute_string(string1) );
		std::string const string2b( EnergyBasedClusteringProtocol::permute_string(string2) );
		std::string const string3b( EnergyBasedClusteringProtocol::permute_string(string3) );
		std::string const string3c( EnergyBasedClusteringProtocol::permute_string(string3b) );
		std::string const string3d( EnergyBasedClusteringProtocol::permute_string(string3c) );
		std::string const string3e( EnergyBasedClusteringProtocol::permute_string(string3d) );

		TS_ASSERT_EQUALS( string1b, "" );
		TS_ASSERT_EQUALS( string2b, "Q" );
		TS_ASSERT_EQUALS( string3b, "orontoT" );
		TS_ASSERT_EQUALS( string3c, "rontoTo" );
		TS_ASSERT_EQUALS( string3d, "ontoTor" );
		TS_ASSERT_EQUALS( string3e, "ntoToro" );
	}

	/// @brief Test that the get_circular_permutation_first_in_alphabetical_order() function works correctly.
	void test_get_circular_permutation_first_in_alphabetical_order() {
		using namespace protocols::energy_based_clustering;

		std::string const input_string1( "XYYXBBAXABYY" );
		std::string const input_string2( "ABYYXYYXBBAX" );
		std::string const input_string3( "ABABABAXABAB" );
		std::string const input_string4( "BABABABAXABA" );

		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::get_circular_permutation_first_in_alphabetical_order(input_string1), "ABYYXYYXBBAX" );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::get_circular_permutation_first_in_alphabetical_order(input_string2), "ABYYXYYXBBAX" );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::get_circular_permutation_first_in_alphabetical_order(input_string3), "ABABABABABAX" );
		TS_ASSERT_EQUALS( EnergyBasedClusteringProtocol::get_circular_permutation_first_in_alphabetical_order(input_string4), "ABABABABABAX" );
	}

	/// @brief Test the EnergyBasedClusteringProtocol::storeposedata() function.
	void test_storeposedata(){
		protocols::energy_based_clustering::EnergyBasedClusteringProtocol clusterer;
		core::pose::Pose pose;
		utility::vector1< core::Real > posedata;
		utility::vector1< numeric::xyzVector <core::Real> > alignmentdata;
		utility::vector1< core::Real > dihedral_mode_resconstruction_data;
		utility::vector1 < core::id::NamedAtomID > extra_atom_list;

		core::pose::make_pose_from_sequence( pose, "AX[OU3_ALA]X[B3A]", "fa_standard", true );
		for ( core::Size i(1); i<=pose.total_residue(); ++i ) {
			pose.set_omega(i, 180);
			pose.set_phi( i, -61 );
			pose.set_psi( i, -41 );
			if ( i > 1 ) pose.set_theta( i, 173.0);
			if ( i == 2 ) pose.set_mu( i, -64.1 );
		}
		pose.update_residue_neighbors();

		clusterer.storeposedata( pose, posedata, alignmentdata, dihedral_mode_resconstruction_data, protocols::energy_based_clustering::EBC_bb_cartesian, extra_atom_list );

		TS_ASSERT_EQUALS( alignmentdata.size(), 17 );

		TR << "\nAlignmentdata:\n";
		for ( core::Size i(1); i<=alignmentdata.size(); ++i ) {
			TR << alignmentdata[i].x() << "\t" << alignmentdata[i].y() << "\t" << alignmentdata[i].z() << "\n";
		}
		TR << std::endl;

		TS_ASSERT_DELTA( alignmentdata[1].x(), pose.xyz( core::id::NamedAtomID( "N", 1 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[1].y(), pose.xyz( core::id::NamedAtomID( "N", 1 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[1].z(), pose.xyz( core::id::NamedAtomID( "N", 1 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[2].x(), pose.xyz( core::id::NamedAtomID( "CA", 1 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[2].y(), pose.xyz( core::id::NamedAtomID( "CA", 1 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[2].z(), pose.xyz( core::id::NamedAtomID( "CA", 1 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[3].x(), pose.xyz( core::id::NamedAtomID( "C", 1 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[3].y(), pose.xyz( core::id::NamedAtomID( "C", 1 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[3].z(), pose.xyz( core::id::NamedAtomID( "C", 1 ) ).z(), 1e-6 );\
			TS_ASSERT_DELTA( alignmentdata[4].x(), pose.xyz( core::id::NamedAtomID( "O", 1 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[4].y(), pose.xyz( core::id::NamedAtomID( "O", 1 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[4].z(), pose.xyz( core::id::NamedAtomID( "O", 1 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[5].x(), pose.xyz( core::id::NamedAtomID( "CB", 1 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[5].y(), pose.xyz( core::id::NamedAtomID( "CB", 1 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[5].z(), pose.xyz( core::id::NamedAtomID( "CB", 1 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[6].x(), pose.xyz( core::id::NamedAtomID( "N", 2 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[6].y(), pose.xyz( core::id::NamedAtomID( "N", 2 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[6].z(), pose.xyz( core::id::NamedAtomID( "N", 2 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[7].x(), pose.xyz( core::id::NamedAtomID( "CA", 2 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[7].y(), pose.xyz( core::id::NamedAtomID( "CA", 2 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[7].z(), pose.xyz( core::id::NamedAtomID( "CA", 2 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[8].x(), pose.xyz( core::id::NamedAtomID( "CM", 2 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[8].y(), pose.xyz( core::id::NamedAtomID( "CM", 2 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[8].z(), pose.xyz( core::id::NamedAtomID( "CM", 2 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[9].x(), pose.xyz( core::id::NamedAtomID( "NU", 2 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[9].y(), pose.xyz( core::id::NamedAtomID( "NU", 2 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[9].z(), pose.xyz( core::id::NamedAtomID( "NU", 2 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[10].x(), pose.xyz( core::id::NamedAtomID( "C", 2 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[10].y(), pose.xyz( core::id::NamedAtomID( "C", 2 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[10].z(), pose.xyz( core::id::NamedAtomID( "C", 2 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[11].x(), pose.xyz( core::id::NamedAtomID( "O", 2 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[11].y(), pose.xyz( core::id::NamedAtomID( "O", 2 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[11].z(), pose.xyz( core::id::NamedAtomID( "O", 2 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[12].x(), pose.xyz( core::id::NamedAtomID( "CB", 2 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[12].y(), pose.xyz( core::id::NamedAtomID( "CB", 2 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[12].z(), pose.xyz( core::id::NamedAtomID( "CB", 2 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[13].x(), pose.xyz( core::id::NamedAtomID( "N", 3 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[13].y(), pose.xyz( core::id::NamedAtomID( "N", 3 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[13].z(), pose.xyz( core::id::NamedAtomID( "N", 3 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[14].x(), pose.xyz( core::id::NamedAtomID( "CA", 3 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[14].y(), pose.xyz( core::id::NamedAtomID( "CA", 3 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[14].z(), pose.xyz( core::id::NamedAtomID( "CA", 3 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[15].x(), pose.xyz( core::id::NamedAtomID( "CM", 3 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[15].y(), pose.xyz( core::id::NamedAtomID( "CM", 3 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[15].z(), pose.xyz( core::id::NamedAtomID( "CM", 3 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[16].x(), pose.xyz( core::id::NamedAtomID( "C", 3 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[16].y(), pose.xyz( core::id::NamedAtomID( "C", 3 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[16].z(), pose.xyz( core::id::NamedAtomID( "C", 3 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[17].x(), pose.xyz( core::id::NamedAtomID( "CB", 3 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[17].y(), pose.xyz( core::id::NamedAtomID( "CB", 3 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[17].z(), pose.xyz( core::id::NamedAtomID( "CB", 3 ) ).z(), 1e-6 );
	}



};
