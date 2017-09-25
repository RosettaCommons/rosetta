// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/rna/AlignmentEnergy.cxxtest.hh
/// @brief  test suite for protocols::rna::AlignmentEnergy
/// @author Andy Watkins (amw579@stanford.edu)

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers

#include <platform/types.hh>

// Package Headers
#include <test/util/deriv_funcs.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/pdb1ubq.hh>
#include <test/core/init_util.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/id/AtomID.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <protocols/rna/AlignmentEnergy.hh>
#include <protocols/stepwise/modeler/align/StepWisePoseAligner.hh>

#include <core/kinematics/FoldTree.hh>
#include <numeric/xyzVector.hh>

//Auto Headers
#include <core/pose/util.hh>
#include <utility/vector1.hh>


// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;

void shift_pose( core::pose::Pose & pose, Real const offset ) {
	using namespace core::id;
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		for ( Size jj = 1; jj <= pose.residue( ii ).natoms(); ++jj ) {
			auto const xyz = pose.xyz( AtomID( jj, ii ) );
			//numeric::xyzVector xyzNew( xyz.x() + offset, xyz.y(), xyz.z() );
			Vector xyzNew( xyz.x() + offset, xyz.y(), xyz.z() );
			pose.set_xyz( AtomID( jj, ii ), xyzNew );

		}
	}
}

class AlignmentEnergyTests : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-rmsd_screen 2.0 -new_align_pdb protocols/rna/srl_fixed_START1_1q9a_RNA.pdb" );
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	/// This tests that the energy works on a partially built pose (not some artificial
	/// shifted finished case)
	void dont_test_VS()
	{
		using namespace core::chemical;
		using namespace core::kinematics;
		using namespace core::pose;
		using namespace core::pose::full_model_info;
		using namespace core::scoring;
		using namespace protocols::rna;

		Pose pose = *core::import_pose::get_pdb_and_cleanup( "protocols/rna/swm_rebuild.out.1.pdb" );
		update_full_model_info_from_pose( pose );

		PoseOP align_pose = core::import_pose::get_pdb_and_cleanup( "protocols/rna/VS_rbzm_P2P3P6_align_ALIGN_4r4v.pdb" );
		update_full_model_info_from_pose( *align_pose );

		core::scoring::methods::EnergyMethodOptions options;
		AlignmentEnergy rna_align( options );
		ScoreFunction sfxn;
		rna_align.align_pose( align_pose );
		rna_align.rmsd_screen( 2.0 );
		EnergyMap totals;
		rna_align.finalize_total_energy( pose, sfxn, totals );
		// This seems to be platform-dependent.
		// 48.5467 or 47.0572.
		//TS_ASSERT_DELTA( totals[ alignment ], 48.5467, 1e-4 );
	}

	/// This can be used to ensure norm matches norm-numeric
	void test_score_values()
	{
		using namespace core::chemical;
		using namespace core::kinematics;
		using namespace core::pose;
		using namespace core::pose::full_model_info;
		using namespace core::scoring;
		using namespace protocols::rna;

		utility::vector1< std::pair< Real, Real > > test_vals;
		test_vals.emplace_back( std::make_pair( 0.0, 0.0 ) );
		test_vals.emplace_back( std::make_pair( 1.0, 0.0 ) );
		test_vals.emplace_back( std::make_pair( 2.0, 0.0 ) );
		test_vals.emplace_back( std::make_pair( 3.0, 1.0 ) );
		test_vals.emplace_back( std::make_pair( 4.0, 4.0 ) );

		for ( Size ii = 1; ii <= test_vals.size(); ++ii ) {

			Real const offset = test_vals[ ii ].first;
			Real const expected = test_vals[ ii ].second;
			Pose pose = *core::import_pose::get_pdb_and_cleanup( "protocols/rna/srl_fixed_START1_1q9a_RNA.pdb" );
			update_full_model_info_from_pose( pose );

			PoseOP align_pose = PoseOP( new Pose( pose ) );
			update_full_model_info_from_pose( *align_pose );

			// Shift all align_pose by one angstrom +x
			shift_pose( pose, offset );

			core::scoring::methods::EnergyMethodOptions options;
			AlignmentEnergy rna_align( options );
			ScoreFunction sfxn;
			rna_align.align_pose( align_pose );
			rna_align.rmsd_screen( 2.0 );
			EnergyMap totals;
			rna_align.finalize_total_energy( pose, sfxn, totals );
			TS_ASSERT_DELTA( totals[ alignment ], expected, 1e-12 );
		}
	}

	void test_deriv() {
		using namespace core::chemical;
		using namespace core::kinematics;
		using namespace core::pose;
		using namespace core::pose::full_model_info;
		using namespace core::scoring;
		using namespace protocols::rna;

		utility::vector1< Real > test_vals;
		test_vals.emplace_back( 3.0 );

		for ( Real const offset : test_vals ) {

			Pose pose = *core::import_pose::get_pdb_and_cleanup( "protocols/rna/srl_fixed_START1_1q9a_RNA.pdb" );
			update_full_model_info_from_pose( pose );

			PoseOP align_pose = PoseOP( new Pose( pose ) );
			update_full_model_info_from_pose( *align_pose );

			// Shift all align_pose by one angstrom +x
			shift_pose( pose, offset );

			ScoreFunction sfxn;
			sfxn.set_weight( fa_atr, 1.0 );
			sfxn.set_weight( fa_rep, 1.0 );
			sfxn.set_weight( fa_sol, 1.0 );
			sfxn.set_weight( fa_elec, 1.0 );
			sfxn.set_weight( alignment, 1.0 );
			sfxn.set_weight( coordinate_constraint, 1.0 );

			MoveMap movemap( create_movemap_to_allow_all_torsions() );
			AtomDerivValidator adv( pose, sfxn, movemap );
			adv.simple_deriv_check( true, 1e-5 );
		}
	}

};


