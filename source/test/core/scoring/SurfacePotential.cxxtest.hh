// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/scoring/SurfacePotential.cxxtest.hh
/// @brief  test suite for core::scoring::SurfacePotential.cc
/// @author Ron Jacak

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/pack/interaction_graph/SurfacePotential.hh>
#include <core/chemical/AA.hh>


#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/TenANeighborGraph.hh>


#include <platform/types.hh>

// Package Headers
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

//Auto Headers
#include <core/id/AtomID.hh>
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>


// --------------- Test Class --------------- //

// using declarations
using namespace core;

class SurfacePotentialTests : public CxxTest::TestSuite {

public:

	pack::interaction_graph::SurfacePotential* sp;
	Real TOLERATED_ERROR;
	pose::Pose pose;
	scoring::ScoreFunctionOP sf;

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
		core_init();
		sp = pack::interaction_graph::SurfacePotential::get_instance();
		TOLERATED_ERROR = 0.01;
	}


	// Shared finalization goes here.
	void tearDown() {
		// sp is just a pointer, and the static class will deallocate when the program ends
	}


	// --------------- Test Cases --------------- //
	void test_average_residue_hASA() {

		// Values from the database file...
		// VAL 97.1757 70.0387 47.4600 23.4852 6.7890
		// TRP 139.783 95.7225 58.8472 34.5932 16.7000
		// ASN 36.3179 27.0235 19.6461 10.9162 4.5142

		TS_ASSERT_DELTA( sp->average_residue_hASA( chemical::aa_val, 3 ), 97.1747, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->average_residue_hASA( chemical::aa_val, 11 ), 70.0387, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->average_residue_hASA( chemical::aa_val, 14 ), 47.4600, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->average_residue_hASA( chemical::aa_val, 18 ), 23.4852, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->average_residue_hASA( chemical::aa_val, 21 ), 6.7890, TOLERATED_ERROR );

		TS_ASSERT_DELTA( sp->average_residue_hASA( chemical::aa_trp, 1 ), 139.783, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->average_residue_hASA( chemical::aa_trp, 12 ), 95.7225, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->average_residue_hASA( chemical::aa_trp, 15 ), 58.8472, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->average_residue_hASA( chemical::aa_trp, 24 ), 16.7000, TOLERATED_ERROR );

		TS_ASSERT_DELTA( sp->average_residue_hASA( chemical::aa_asn, 10 ), 36.3179, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->average_residue_hASA( chemical::aa_asn, 13 ), 27.0235, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->average_residue_hASA( chemical::aa_asn, 16 ), 19.6461, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->average_residue_hASA( chemical::aa_asn, 19 ), 10.9162, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->average_residue_hASA( chemical::aa_asn, 23 ), 4.5142, TOLERATED_ERROR );
	}


	/// @brief old version of this test that uses a cutoff of 16nbs only for surface exposed residues
	void x_test_hASA_patch_energy() {

		// Values from the database file... 0-24 is 3.524; 25-50 is something else...
		// 0  0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
		// 200  0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
		// 500  2.785 2.785 2.785 2.785 2.785 2.378 1.883 1.340 1.160 0.732 0.617 0.415 0.118 0.184 0.039 0.000
		// 775  25.000 25.000 25.000 25.000 25.000 25.000 25.000 6.402 5.203 4.464 4.467 4.549 4.057 3.812 3.309 2.766
		// 900  25.000 25.000 25.000 25.000 25.000 25.000 25.000 25.000 6.812 7.029 6.546 5.530 5.561 5.758 5.137 5.091

		TS_ASSERT_DELTA( sp->hASA_patch_energy( 5, 1 ), 0.000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->hASA_patch_energy( 212, 5 ), 0.000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->hASA_patch_energy( 501, 1 ), 2.785, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->hASA_patch_energy( 501, 10 ), 0.732, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->hASA_patch_energy( 775, 5 ), 25.000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->hASA_patch_energy( 775, 15 ), 3.309, TOLERATED_ERROR );
		//TS_ASSERT_DELTA( sp->hASA_patch_energy( 900, 18 ), 5.758, TOLERATED_ERROR );
	}

	void test_hASA_patch_energy() {

		// Values from the database file... 0-24 is 3.524; 25-50 is something else...
		// 0 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
		// 200 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
		// 500 2.785 2.785 2.785 2.785 2.785 2.378 1.883 1.340 1.160 0.732 0.617 0.415 0.118 0.184 0.039 0.000 0.000 0.008 0.000 0.053
		// 775 25.000 25.000 25.000 25.000 25.000 25.000 25.000 6.402 5.203 4.464 4.467 4.549 4.057 3.812 3.309 2.766 2.687 2.175 1.953 1.780
		// 900 25.000 25.000 25.000 25.000 25.000 25.000 25.000 25.000 6.812 7.029 6.546 5.530 5.561 5.758 5.137 5.091 4.938 4.788 4.027 3.787

		TS_ASSERT_DELTA( sp->hASA_patch_energy( 5, 1 ), 0.000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->hASA_patch_energy( 212, 5 ), 0.000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->hASA_patch_energy( 501, 1 ), 2.785, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->hASA_patch_energy( 501, 10 ), 0.732, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->hASA_patch_energy( 501, 20 ), 0.053, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->hASA_patch_energy( 775, 5 ), 25.000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->hASA_patch_energy( 775, 15 ), 3.309, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->hASA_patch_energy( 775, 20 ), 1.780, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->hASA_patch_energy( 900, 18 ), 4.788, TOLERATED_ERROR );
	}

	void test_compute_residue_surface_energy() {

		core::import_pose::pose_from_pdb( pose, "core/pack/1l2y_renameH.pdb" );

		// create a score function using the standard packer weights
		sf = scoring::get_score_function();
		(*sf)( pose );

		scoring::EnergyMap emap;

		utility::vector1< Size > nbs;
		nbs.resize( pose.n_residue(), 0 );
		scoring::TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );
		for ( Size ii=1; ii <= pose.n_residue(); ++ii ) {
			nbs[ ii ] = tenA_neighbor_graph.get_node( ii )->num_neighbors_counting_self();
		}

		emap.zero();
		sp->compute_residue_surface_energy( pose.residue(6), pose, emap, 6, nbs);
		TS_ASSERT_DELTA( emap[ scoring::surface ], 25.0000, TOLERATED_ERROR );

		emap.zero();
		sp->compute_residue_surface_energy( pose.residue(7), pose, emap, 7, nbs);
		TS_ASSERT_DELTA( emap[ scoring::surface ], 3.6329, TOLERATED_ERROR );

		emap.zero();
		sp->compute_residue_surface_energy( pose.residue(8), pose, emap, 8, nbs);
		TS_ASSERT_DELTA( emap[ scoring::surface ], 1.7410, TOLERATED_ERROR );
	}

	void test_compute_pose_surface_energy_fragment() {

		chemical::ResidueTypeSetCOP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
		core::pose::make_pose_from_sequence( pose, "DFGLKANM", *rsd_set );

		for ( Size ii=1; ii <= pose.n_residue(); ii+=3 ) {
			pose.set_phi( ii, -150.0 );
			pose.set_psi( ii, 150.0 );
			pose.set_omega( ii, 180.0 );
		}

		Real total_surfaceE = -1.0;
		utility::vector1< Real > residue_surfaceE( pose.n_residue(), 0.0 );

		sp->compute_pose_surface_energy( pose, total_surfaceE, residue_surfaceE );
		TS_ASSERT_DELTA( total_surfaceE, 0.0000, TOLERATED_ERROR );  // tests unscored poses; return 0.0 for all residues

		sf = scoring::get_score_function();
		(*sf)( pose );

		sp->compute_pose_surface_energy( pose, total_surfaceE, residue_surfaceE );
		TS_ASSERT_DELTA( total_surfaceE, 10.2460, TOLERATED_ERROR );

		TS_ASSERT_DELTA( residue_surfaceE[1], 2.3780, TOLERATED_ERROR );
		TS_ASSERT_DELTA( residue_surfaceE[4], 1.9600, TOLERATED_ERROR );
		TS_ASSERT_DELTA( residue_surfaceE[5], 1.9600, TOLERATED_ERROR );
		TS_ASSERT_DELTA( residue_surfaceE[6], 0.3790, TOLERATED_ERROR );
	}

	void test_compute_pose_surface_energy_smallprotein() {

		core::import_pose::pose_from_pdb( pose, "core/pack/1l2y_renameH.pdb" );

		Real total_surfaceE = -1.0;
		utility::vector1< Real > residue_surfaceE( pose.n_residue(), 0.0 );

		sp->compute_pose_surface_energy( pose, total_surfaceE, residue_surfaceE );
		TS_ASSERT_DELTA( total_surfaceE, 0.0000, TOLERATED_ERROR );  // tests unscored poses; return 0.0 for all residues

		sf = scoring::get_score_function();
		(*sf)( pose );

		sp->compute_pose_surface_energy( pose, total_surfaceE, residue_surfaceE );
		TS_ASSERT_DELTA( total_surfaceE, 67.8569, TOLERATED_ERROR );

		TS_ASSERT_DELTA( residue_surfaceE[1], 1.8829, TOLERATED_ERROR );
		TS_ASSERT_DELTA( residue_surfaceE[4], 3.8509, TOLERATED_ERROR );
		TS_ASSERT_DELTA( residue_surfaceE[5], 3.8509, TOLERATED_ERROR );
		TS_ASSERT_DELTA( residue_surfaceE[6], 25.0000, TOLERATED_ERROR );
	}

	void test_hpatch_score() {

		// Values from the database file...
		// 0 0.000
		// 50 0.160
		// 400 10.240
		// 850 46.240
		// 900 100.000

		TS_ASSERT_DELTA( sp->hpatch_score( 5 ), 0.0, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->hpatch_score( 51 ), 0.160, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->hpatch_score( 401 ), 10.240, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->hpatch_score( 860 ), 46.240, TOLERATED_ERROR );
		TS_ASSERT_DELTA( sp->hpatch_score( 924 ), 100.00, TOLERATED_ERROR );
	}

	void test_compute_pose_hpatch_energy() {

		core::import_pose::pose_from_pdb( pose, "core/pack/1l2y_renameH.pdb" );

		Real total_hpatch_score = 0.0;
		std::map< Size, std::pair< Real, Real > > patch_scores;
		std::map< Size, utility::vector1< id::AtomID > > atoms_in_patches;
		sp->compute_pose_hpatch_score( pose, total_hpatch_score, patch_scores, atoms_in_patches );

		TS_ASSERT_DELTA( total_hpatch_score, 19.3600, TOLERATED_ERROR );
		TS_ASSERT_DELTA( patch_scores[ 16 ].first, 19.3600, TOLERATED_ERROR );

		// test a more realistic pose than Trp-cage
		pose = create_1ten_pdb_pose();

		total_hpatch_score = 0.0;
		patch_scores.clear();
		atoms_in_patches.clear();
		sp->compute_pose_hpatch_score( pose, total_hpatch_score, patch_scores, atoms_in_patches );

		TS_ASSERT_DELTA( total_hpatch_score, 4.9600, TOLERATED_ERROR );
		TS_ASSERT_DELTA( patch_scores[ 170 ].first, 0.160, TOLERATED_ERROR );
		TS_ASSERT_DELTA( patch_scores[ 196 ].first, 0.6400, TOLERATED_ERROR );
		TS_ASSERT_DELTA( patch_scores[ 289 ].first, 4.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( patch_scores[ 617 ].first, 0.1600, TOLERATED_ERROR );
		TS_ASSERT_DELTA( patch_scores[ 680 ].first, 0.00, TOLERATED_ERROR );

	}

	void test_compute_pose_with_HETATM_hpatch_energy() {

		// the following PDB has Ca and some other HETATMs in it
		core::import_pose::pose_from_pdb( pose, "core/pack/2mcm_0001.pdb" );

		Real total_hpatch_score = 0.0;
		std::map< Size, std::pair< Real, Real > > patch_scores;
		std::map< Size, utility::vector1< id::AtomID > > atoms_in_patches;
		sp->compute_pose_hpatch_score( pose, total_hpatch_score, patch_scores, atoms_in_patches );

		// checks that the function ignores HETATMs correctly, and comes up with the right score
		TS_ASSERT_DELTA( total_hpatch_score, 29.6000, TOLERATED_ERROR );
	}

};

