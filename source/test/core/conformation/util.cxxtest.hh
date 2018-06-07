// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    test/core/chemical/rings/util.cxxtest.hh
/// @brief   Test suite for conformation-related utility functions
/// @author  Labonte <JWLabonte@jhu.edu>


// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Header
#include <core/conformation/util.hh>

// Package Header
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/carbohydrates/util.hh>

// Project Headers
#include <core/types.hh>
#include <core/chemical/rings/RingConformer.hh>
#include <core/chemical/rings/RingConformerSet.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/id/TorsionID.hh>

// Utility Header
#include <utility/vector1.hh>

// Basic Header
#include <basic/Tracer.hh>


static basic::Tracer TR( "core.conformation.util.cxxtest.hh" );


class ConformationUtilityFunctionTests : public CxxTest::TestSuite {
public:  // Standard methods //////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		core_init_with_additional_options( "-include_sugars" );
	}

	// Destruction
	void tearDown()
	{}


public:  // Tests /////////////////////////////////////////////////////////////
	void test_find_bond_torsion_with_nearest_orientation()
	{
		using namespace std;
		using namespace utility;
		using namespace core::id;
		using namespace core::chemical::rings;
		using namespace core::pose;
		using namespace core::conformation;
		using namespace core::conformation::carbohydrates;

		TR << "Testing that the best torsions are chosen to minimize downstream effects of a move." << endl;
		TR << "Testing a peptide." << endl;

		Pose peptide;
		make_pose_from_sequence( peptide, "AAA", "fa_standard" );
		// Ensure that these are trans peptides.
		peptide.set_omega( 1, 180.0 );
		peptide.set_omega( 2, 180.0 );

		vector1< TorsionID > const peptide_torsions = {
			TorsionID( 1, BB, 1 ),  // Phi1
			TorsionID( 1, BB, 2 ),  // Psi1
			TorsionID( 2, BB, 1 ),  // Phi2
			TorsionID( 2, BB, 2 ),  // Psi2
			TorsionID( 3, BB, 1 ),  // Phi3
			TorsionID( 3, BB, 2 )   // Psi3
			};

		TS_ASSERT_EQUALS( find_bond_torsion_with_nearest_orientation(
			peptide.conformation(), peptide_torsions, TorsionID( 2, BB, 1 ) ),
			TorsionID( 1, BB, 2 ) );  // Because omega is trans, Psi1 should undo a move at Phi2.
		TS_ASSERT_EQUALS( find_bond_torsion_with_nearest_orientation(
			peptide.conformation(), peptide_torsions, TorsionID( 2, BB, 2 ) ),
			TorsionID( 3, BB, 1 ) );  // Because omega is trans, Phi3 should undo a move at Psi2.


		vector1< TorsionID > sugar_torsions;
		TorsionID phi2;
		TorsionID psi2;
		TorsionID omega2;
		TorsionID phi3;
		TorsionID psi3;
		TorsionID omega3;
		TorsionID phi4;
		TorsionID psi4;
		TorsionID omega4;


		TR << "Testing a beta-1,4-oligosaccharide with 4C1 all-equatorial rings." << endl;

		Pose to4_beta_sugar;
		make_pose_from_saccharide_sequence( to4_beta_sugar,
			"b-D-Glcp-(1->4)-b--D-Glcp-(1->4)-b-D-Glcp-(1->4)-b-D-Glcp", "fa_standard" );

		phi2 = get_non_NU_TorsionID_from_AtomIDs(
			to4_beta_sugar.conformation(), get_reference_atoms_for_phi( to4_beta_sugar.conformation(), 2 ) );
		psi2 = get_non_NU_TorsionID_from_AtomIDs(
			to4_beta_sugar.conformation(), get_reference_atoms_for_psi( to4_beta_sugar.conformation(), 2 ) );
		phi3 = get_non_NU_TorsionID_from_AtomIDs(
			to4_beta_sugar.conformation(), get_reference_atoms_for_phi( to4_beta_sugar.conformation(), 3 ) );
		psi3 = get_non_NU_TorsionID_from_AtomIDs(
			to4_beta_sugar.conformation(), get_reference_atoms_for_psi( to4_beta_sugar.conformation(), 3 ) );
		phi4 = get_non_NU_TorsionID_from_AtomIDs(
			to4_beta_sugar.conformation(), get_reference_atoms_for_phi( to4_beta_sugar.conformation(), 4 ) );
		psi4 = get_non_NU_TorsionID_from_AtomIDs(
			to4_beta_sugar.conformation(), get_reference_atoms_for_psi( to4_beta_sugar.conformation(), 4 ) );

		sugar_torsions.push_back( phi2 );
		sugar_torsions.push_back( psi2 );
		sugar_torsions.push_back( phi3 );
		sugar_torsions.push_back( psi3 );
		sugar_torsions.push_back( phi4 );
		sugar_torsions.push_back( psi4 );

		TS_ASSERT_EQUALS( find_bond_torsion_with_nearest_orientation(
			to4_beta_sugar.conformation(), sugar_torsions, phi3),
			psi4 );  // 4C1 beta-glucose is all-equatorial; the ring holds Psi4 near-parallel to Phi3.


		TR << "Testing a alpha-1,6-oligosaccharide with 1C4 all-axial rings." << endl;

		Pose to6_beta_sugar;
		make_pose_from_saccharide_sequence( to6_beta_sugar,
			"b-D-Glcp-(1->6)-b--D-Glcp-(1->6)-b-D-Glcp-(1->6)-b-D-Glcp", "fa_standard" );
		RingConformerSet const conformers( 6, "", vector1< string >() );
		RingConformer const & flipped_chair( conformers.get_ideal_conformer_by_name( "1C4" ) );

		for ( core::uint i( 1 ); i <= 4; ++i ) {
			to6_beta_sugar.set_psi( i, 90.0 );  // Ensure that none of the psis are near 180 so that test passes.
			to6_beta_sugar.set_ring_conformation( i, 1, flipped_chair );
		}

		phi2 = get_non_NU_TorsionID_from_AtomIDs(
			to6_beta_sugar.conformation(), get_reference_atoms_for_phi( to6_beta_sugar.conformation(), 2 ) );
		psi2 = get_non_NU_TorsionID_from_AtomIDs(
			to6_beta_sugar.conformation(), get_reference_atoms_for_psi( to6_beta_sugar.conformation(), 2 ) );
		omega2 = get_non_NU_TorsionID_from_AtomIDs(
			to6_beta_sugar.conformation(), get_reference_atoms_for_1st_omega( to6_beta_sugar.conformation(), 2 ) );
		phi3 = get_non_NU_TorsionID_from_AtomIDs(
			to6_beta_sugar.conformation(), get_reference_atoms_for_phi( to6_beta_sugar.conformation(), 3 ) );
		psi3 = get_non_NU_TorsionID_from_AtomIDs(
			to6_beta_sugar.conformation(), get_reference_atoms_for_psi( to6_beta_sugar.conformation(), 3 ) );
		omega3 = get_non_NU_TorsionID_from_AtomIDs(
			to6_beta_sugar.conformation(), get_reference_atoms_for_1st_omega( to6_beta_sugar.conformation(), 3 ) );
		phi4 = get_non_NU_TorsionID_from_AtomIDs(
			to6_beta_sugar.conformation(), get_reference_atoms_for_phi( to6_beta_sugar.conformation(), 4 ) );
		psi4 = get_non_NU_TorsionID_from_AtomIDs(
			to6_beta_sugar.conformation(), get_reference_atoms_for_psi( to6_beta_sugar.conformation(), 4 ) );
		omega4 = get_non_NU_TorsionID_from_AtomIDs(
			to6_beta_sugar.conformation(), get_reference_atoms_for_1st_omega( to6_beta_sugar.conformation(), 4 ) );

		sugar_torsions.clear();
		sugar_torsions.push_back( phi2 );
		sugar_torsions.push_back( psi2 );
		sugar_torsions.push_back( omega2 );
		sugar_torsions.push_back( phi3 );
		sugar_torsions.push_back( psi3 );
		sugar_torsions.push_back( omega3 );
		sugar_torsions.push_back( phi4 );
		sugar_torsions.push_back( psi4 );
		sugar_torsions.push_back( omega4 );

		TS_ASSERT_EQUALS( find_bond_torsion_with_nearest_orientation(
			to6_beta_sugar.conformation(), sugar_torsions, phi3 ),
			omega4 );  // 1C4 beta-glucose is all-axial; the ring holds Omega4 near-parallel to Phi3.
	}

	// Confirm that the position of a substituent on a ring can be correctly assigned.
	void test_position_of_atom_on_ring()
	{
		using namespace core::conformation;
		using namespace core::pose;

		TR << "Testing that that position_of_atom_on_ring() correctly figures out attachment positions." << std::endl;

		Pose monosaccharide;
		make_pose_from_saccharide_sequence( monosaccharide, "alpha-D-Glcp", "fa_standard" );
		Residue const & glucose( monosaccharide.residue( 1 ) );
		utility::vector1< core::uint > const & ring_atoms( glucose.type().ring_atoms( 1 ) );

		uint const O1_index( glucose.atom_index( " O1 " ) );
		uint const O2_index( glucose.atom_index( " O2 " ) );
		uint const O3_index( glucose.atom_index( " O3 " ) );
		uint const O4_index( glucose.atom_index( " O4 " ) );
		uint const O5_index( glucose.atom_index( " O5 " ) );
		uint const O6_index( glucose.atom_index( " O6 " ) );

		TS_ASSERT_EQUALS( position_of_atom_on_ring( glucose, O1_index, ring_atoms ), 1 );
		TS_ASSERT_EQUALS( position_of_atom_on_ring( glucose, O2_index, ring_atoms ), 2 );
		TS_ASSERT_EQUALS( position_of_atom_on_ring( glucose, O3_index, ring_atoms ), 3 );
		TS_ASSERT_EQUALS( position_of_atom_on_ring( glucose, O4_index, ring_atoms ), 4 );
		TS_ASSERT_EQUALS( position_of_atom_on_ring( glucose, O5_index, ring_atoms ), 5 );  // O5 is in the ring itself.
		TS_ASSERT_EQUALS( position_of_atom_on_ring( glucose, O6_index, ring_atoms ), 0 );  // O6 is exocyclic.
	}


	void test_residue_from_name() {
		using namespace core::conformation;
		ResidueCOP res = get_residue_from_name1( 'Q', true, false, false );

		TS_ASSERT( res->name1() == 'Q' );
		TS_ASSERT( res->type().has_variant_type( core::chemical::LOWER_TERMINUS_VARIANT ) );

		std::string name = res->name();
		res = get_residue_from_name( name );

		TS_ASSERT( res->name1() == 'Q' );
		TS_ASSERT( res->type().has_variant_type( core::chemical::LOWER_TERMINUS_VARIANT ) );

		utility::vector1< core::chemical::VariantType > variant { core::chemical::REPLONLY };
		res = get_residue_from_name1( 'W', false, true, true );

		TS_ASSERT( res->name1() == 'W' );
		TS_ASSERT( res->type().has_variant_type( core::chemical::UPPER_TERMINUS_VARIANT ) );
		TS_ASSERT( res->type().is_d_aa() );
		TR << "Residue from name passed" << std::endl;
	}

};  // class ConformationUtilityFunctionTests
