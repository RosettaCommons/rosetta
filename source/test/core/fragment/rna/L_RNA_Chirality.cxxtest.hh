// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/fragment/rna/L_RNA_Chirality.cxxtest.hh
/// @brief  Mirror-image (L-RNA) handling: fragment insertion, pucker helpers, ideal geometry, scoring.
/// @details Each test builds a D-RNA pose and its exact mirror image out of the built-in L-RNA
///  types (one-letter codes 1-4), then checks that L-RNA code paths treat the mirror image exactly
///  as the D-RNA paths treat the original.

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/rna/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/fragment/rna/FullAtomRNA_Fragments.hh>
#include <core/fragment/rna/TorsionSet.hh>
#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/rna/RNA_IdealCoord.hh>
#include <core/pose/toolbox/AtomLevelDomainMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>

#include <basic/database/open.hh>
#include <basic/Tracer.hh>

#include <cmath>

static basic::Tracer TR( "core.fragment.rna.L_RNA_Chirality.cxxtest" );

using namespace core;
using namespace core::chemical::rna;

class L_RNA_ChiralityTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	// All RNA torsions (backbone + chi, nu2, nu1, 2'-OH) that exist in this pose.
	utility::vector1< id::TorsionID >
	rna_torsion_ids( pose::Pose const & pose ) {
		utility::vector1< id::TorsionID > torsion_ids;
		for ( Size i = 1; i <= pose.size(); ++i ) {
			for ( Size j = 1; j <= NUM_RNA_TORSIONS; ++j ) {
				id::TorsionID const tid = ( j <= NUM_RNA_MAINCHAIN_TORSIONS ) ?
					id::TorsionID( i, id::BB, j ) : id::TorsionID( i, id::CHI, j - NUM_RNA_MAINCHAIN_TORSIONS );
				id::AtomID id1, id2, id3, id4;
				if ( pose.conformation().get_torsion_angle_atom_ids( tid, id1, id2, id3, id4 ) ) continue; // true means failure
				torsion_ids.push_back( tid );
			}
		}
		return torsion_ids;
	}

	// Make l_pose the mirror image of d_pose by negating every torsion.
	void
	mirror_torsions( pose::Pose const & d_pose, pose::Pose & l_pose ) {
		for ( id::TorsionID const & tid : rna_torsion_ids( d_pose ) ) {
			l_pose.set_torsion( tid, -d_pose.torsion( tid ) );
		}
	}

	// Mirror images have identical internal distances.
	void
	assert_mirror_images( pose::Pose const & d_pose, pose::Pose const & l_pose, Real const tol = 1.0e-3 ) {
		TS_ASSERT_EQUALS( d_pose.size(), l_pose.size() );
		for ( Size i = 1; i <= d_pose.size(); ++i ) {
			TS_ASSERT_EQUALS( d_pose.residue( i ).natoms(), l_pose.residue( i ).natoms() );
			for ( Size j = i; j <= d_pose.size(); ++j ) {
				for ( Size ii = 1; ii <= d_pose.residue( i ).nheavyatoms(); ++ii ) {
					for ( Size jj = 1; jj <= d_pose.residue( j ).nheavyatoms(); ++jj ) {
						Real const d_dist = d_pose.residue( i ).xyz( ii ).distance( d_pose.residue( j ).xyz( jj ) );
						Real const l_dist = l_pose.residue( i ).xyz( ii ).distance( l_pose.residue( j ).xyz( jj ) );
						TS_ASSERT_DELTA( d_dist, l_dist, tol );
					}
				}
			}
		}
	}

	void test_one_letter_codes_give_l_rna() {
		pose::Pose l_pose;
		pose::make_pose_from_sequence( l_pose, "1234", "fa_standard" );
		TS_ASSERT_EQUALS( l_pose.sequence(), "1234" );
		for ( Size i = 1; i <= l_pose.size(); ++i ) TS_ASSERT( l_pose.residue_type( i ).is_l_rna() );
		TS_ASSERT_EQUALS( l_pose.residue_type( 1 ).aa(), chemical::na_lra );
		TS_ASSERT_EQUALS( l_pose.residue_type( 2 ).aa(), chemical::na_lrc );
		TS_ASSERT_EQUALS( l_pose.residue_type( 3 ).aa(), chemical::na_lrg );
		TS_ASSERT_EQUALS( l_pose.residue_type( 4 ).aa(), chemical::na_lur );
		TS_ASSERT_EQUALS( chemical::oneletter_code_from_aa( chemical::na_lrg ), '3' );
	}

	void test_mirrored_types_and_torsions_give_mirror_image() {
		pose::Pose d_pose, l_pose;
		pose::make_pose_from_sequence( d_pose, "gcaugc", "fa_standard" );
		pose::make_pose_from_sequence( l_pose, "321432", "fa_standard" );
		pose::rna::apply_ideal_c2endo_sugar_coords( d_pose, 3 );
		mirror_torsions( d_pose, l_pose );
		assert_mirror_images( d_pose, l_pose );
	}

	void test_pucker_assignment_and_ideal_c2endo() {
		pose::Pose d_pose, l_pose;
		pose::make_pose_from_sequence( d_pose, "gcaugc", "fa_standard" );
		pose::make_pose_from_sequence( l_pose, "321432", "fa_standard" );
		pose::rna::apply_ideal_c2endo_sugar_coords( d_pose, 3 );
		pose::rna::apply_ideal_c2endo_sugar_coords( l_pose, 3 );
		TS_ASSERT_DELTA( l_pose.torsion( id::TorsionID( 3, id::BB, DELTA ) ), -d_pose.torsion( id::TorsionID( 3, id::BB, DELTA ) ), 1.0e-3 );
		for ( Size i = 1; i <= d_pose.size(); ++i ) {
			TS_ASSERT_EQUALS( pose::rna::assign_pucker( d_pose, i ), pose::rna::assign_pucker( l_pose, i ) );
			TS_ASSERT_EQUALS( get_residue_pucker_state( d_pose.residue( i ) ), get_residue_pucker_state( l_pose.residue( i ) ) );
			TS_ASSERT_EQUALS( get_residue_base_state( d_pose.residue( i ) ), get_residue_base_state( l_pose.residue( i ) ) );
		}
		TS_ASSERT_EQUALS( pose::rna::assign_pucker( l_pose, 3 ), SOUTH );
		// The rest of each pose is still the default A-form, so the two are full mirror images.
		mirror_torsions( d_pose, l_pose );
		assert_mirror_images( d_pose, l_pose );
	}

	void test_ideal_coord_pucker_is_mirrored() {
		pose::Pose d_pose, l_pose;
		pose::make_pose_from_sequence( d_pose, "gcaugc", "fa_standard" );
		pose::make_pose_from_sequence( l_pose, "321432", "fa_standard" );
		mirror_torsions( d_pose, l_pose );
		pose::rna::RNA_IdealCoord const ideal_coord;
		for ( Size i = 1; i <= d_pose.size(); ++i ) {
			ideal_coord.apply( d_pose, i, SOUTH );
			ideal_coord.apply( l_pose, i, SOUTH );
		}
		for ( Size i = 1; i <= d_pose.size(); ++i ) {
			TS_ASSERT_EQUALS( pose::rna::assign_pucker( d_pose, i ), SOUTH );
			TS_ASSERT_EQUALS( pose::rna::assign_pucker( l_pose, i ), SOUTH );
		}
		assert_mirror_images( d_pose, l_pose );
	}

	void test_fragment_insertion_inverts_for_l_and_skips_mixed_windows() {
		using namespace core::fragment::rna;
		FullAtomRNA_Fragments const fragments( basic::database::full_name( "sampling/rna/1jj2.torsions" ) );

		// A two-residue fragment with a C2'-endo first residue and non-A-form backbone.
		TorsionSet torsion_set( 2, 1 );
		Real const values[ 2 ][ NUM_RNA_TORSIONS ] = {
			{ 69.4, -173.0, 58.9, 147.2, -85.4, -38.4, 111.7, -36.4, 156.4, 179.9 },
			{ -68.0, 178.0, 54.0, 82.0, -150.0, -71.0, 79.0, 37.0, -25.0, 100.0 } };
		for ( Size offset = 0; offset < 2; ++offset ) {
			for ( Size j = 1; j <= NUM_RNA_TORSIONS; ++j ) torsion_set.torsions( j, offset ) = values[ offset ][ j - 1 ];
		}

		pose::Pose d_pose, l_pose, mixed_pose;
		pose::make_pose_from_sequence( d_pose, "gcaugc", "fa_standard" );
		pose::make_pose_from_sequence( l_pose, "321432", "fa_standard" );
		pose::make_pose_from_sequence( mixed_pose, "gca432", "fa_standard" );

		pose::toolbox::AtomLevelDomainMapOP d_map( new pose::toolbox::AtomLevelDomainMap( d_pose ) );
		pose::toolbox::AtomLevelDomainMapOP l_map( new pose::toolbox::AtomLevelDomainMap( l_pose ) );
		pose::toolbox::AtomLevelDomainMapOP mixed_map( new pose::toolbox::AtomLevelDomainMap( mixed_pose ) );
		d_map->set( true ); l_map->set( true ); mixed_map->set( true );

		bool is_l( false );
		TS_ASSERT( window_has_uniform_rna_chirality( d_pose, 3, 2, is_l ) && !is_l );
		TS_ASSERT( window_has_uniform_rna_chirality( l_pose, 3, 2, is_l ) && is_l );
		TS_ASSERT( !window_has_uniform_rna_chirality( mixed_pose, 3, 2, is_l ) );

		fragments.insert_fragment( d_pose, 3, torsion_set, d_map );
		fragments.insert_fragment( l_pose, 3, torsion_set, l_map );
		for ( Size offset = 0; offset < 2; ++offset ) {
			for ( Size j = 1; j <= NUM_RNA_MAINCHAIN_TORSIONS; ++j ) {
				id::TorsionID const tid( 3 + offset, id::BB, j );
				TS_ASSERT_DELTA( d_pose.torsion( tid ), values[ offset ][ j - 1 ], 1.0e-3 );
				TS_ASSERT_DELTA( l_pose.torsion( tid ), -values[ offset ][ j - 1 ], 1.0e-3 );
			}
			id::TorsionID const chi( 3 + offset, id::CHI, 1 );
			TS_ASSERT_DELTA( l_pose.torsion( chi ), -d_pose.torsion( chi ), 1.0e-3 );
		}
		mirror_torsions( d_pose, l_pose ); // only touches torsions; should already hold
		assert_mirror_images( d_pose, l_pose );

		// Mixed window: nothing changes.
		pose::Pose const mixed_before( mixed_pose );
		fragments.insert_fragment( mixed_pose, 3, torsion_set, mixed_map );
		for ( id::TorsionID const & tid : rna_torsion_ids( mixed_pose ) ) {
			TS_ASSERT_DELTA( mixed_pose.torsion( tid ), mixed_before.torsion( tid ), 1.0e-6 );
		}
	}

	// Scoring a mirror-image L-RNA pose must give the same energies, term by term, as the D pose.
	void
	check_mirror_energies( std::string const & weights ) {
		scoring::ScoreFunctionOP sfxn = scoring::ScoreFunctionFactory::create_score_function( weights );
		sfxn->set_weight( scoring::rna_suite, 1.0 );
		sfxn->set_weight( scoring::rna_jr_suite, 1.0 );
		sfxn->set_weight( scoring::rna_sugar_close, 1.0 );
		sfxn->set_weight( scoring::rna_torsion, 1.0 );

		pose::Pose d_pose, l_pose;
		pose::make_pose_from_sequence( d_pose, "ggcuucggcc", "fa_standard" );
		pose::make_pose_from_sequence( l_pose, "3324433422", "fa_standard" );
		// Some non-A-form structure so that suites, puckers and stacking are not all identical.
		// (Applied to both poses: the C2'-endo sugar ring has non-ideal bond geometry, not just torsions.)
		for ( Size const i : { 5, 6 } ) {
			pose::rna::apply_ideal_c2endo_sugar_coords( d_pose, i );
			pose::rna::apply_ideal_c2endo_sugar_coords( l_pose, i );
		}
		d_pose.set_torsion( id::TorsionID( 6, id::BB, ALPHA ), 70.0 );
		d_pose.set_torsion( id::TorsionID( 7, id::CHI, 1 ), -60.0 ); // syn
		mirror_torsions( d_pose, l_pose );
		assert_mirror_images( d_pose, l_pose );

		Real const d_total = ( *sfxn )( d_pose );
		Real const l_total = ( *sfxn )( l_pose );
		for ( scoring::ScoreType const st : sfxn->get_nonzero_weighted_scoretypes() ) {
			Real const d_e = d_pose.energies().total_energies()[ st ];
			Real const l_e = l_pose.energies().total_energies()[ st ];
			TR << weights << " " << scoring::name_from_score_type( st ) << " D: " << d_e << " L: " << l_e << std::endl;
			TS_ASSERT_DELTA( l_e, d_e, 1.0e-3 * std::max( 1.0, std::abs( d_e ) ) );
		}
		TS_ASSERT_DELTA( l_total, d_total, 1.0e-3 * std::max( 1.0, std::abs( d_total ) ) );
	}

	void test_mirror_energies_hires() {
		check_mirror_energies( "stepwise/rna/rna_res_level_energy4" );
	}

	void test_mirror_energies_lores() {
		check_mirror_energies( "rna/denovo/rna_lores" );
	}

};
