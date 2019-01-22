// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/lkball/LK_BallEnergy.cxxtest.hh
/// @brief  Unit tests for the lk-ball energy function
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/deriv_funcs.hh>

#include <core/types.hh>

// Unit headers
#include <core/scoring/lkball/LK_BallEnergy.hh>
#include <core/scoring/lkball/LK_BallInfo.hh>

// Package headers
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/trie/RotamerTrie.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/MinimizationData.hh>

// Project headers
#include <core/conformation/residue_datacache.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

#include <utility/graph/Graph.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/conformation/AbstractRotamerTrie.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/backtrace.hh>

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION



using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.lkball.LK_BallEnergy.cxxtest");

// using declarations
using namespace std;
using namespace core;
using namespace scoring;
using namespace etable;
using namespace lkball;


///////////////////////////////////////////////////////////////////////////
/// @name EtableEnergyTest
/// @brief: Test the functionality of the EtableEnergy class
///////////////////////////////////////////////////////////////////////////
class LK_BallEnergyTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	LKB_ResidueInfo const &
	retrieve_lkb_resdata(
		conformation::Residue const & res
	)
	{
		using namespace core::conformation::residue_datacache;
		debug_assert( utility::pointer::dynamic_pointer_cast< LKB_ResidueInfo const > ( res.data().get_const_ptr( LK_BALL_INFO )));
		return ( static_cast< LKB_ResidueInfo const & > ( res.data().get( LK_BALL_INFO ) ) );
	}

	void test_lk_ball_setup_for_scoring() {
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( core::scoring::lk_ball, 1.0 );
		sfxn( pose );

		core::conformation::Residue const & res1 = pose.residue(1);
		LKB_ResidueInfo const & r1_lkbinfo = retrieve_lkb_resdata(res1);

		TS_ASSERT_EQUALS( res1.natoms(), r1_lkbinfo.n_attached_waters().size() );
		TS_ASSERT_EQUALS( res1.natoms(), r1_lkbinfo.water_offset_for_atom().size() );

		utility::vector1< Size > actual_offsets( res1.natoms(), 0 );
		for ( Size ii = 1; ii < res1.natoms(); ++ii ) {
			actual_offsets[ ii+1 ] = actual_offsets[ ii ] + r1_lkbinfo.n_attached_waters(ii);
		}
		Size n_waters_actual = actual_offsets[ res1.natoms() ] + r1_lkbinfo.n_attached_waters(res1.natoms());
		TS_ASSERT_EQUALS( r1_lkbinfo.waters().size(), n_waters_actual );
		for ( Size ii = 1; ii <= res1.nheavyatoms(); ++ii ) {
			TS_ASSERT_EQUALS( actual_offsets[ ii ], r1_lkbinfo.water_offset_for_atom(ii) );
			TS_ASSERT_EQUALS( actual_offsets[ ii ], r1_lkbinfo.water_offset_for_atom()[ii] );
		}

#ifndef NDEBUG
		// The derivative matrices are not allocated until setup_for_derivatives is called
		// so we should get an assertion failure here
		// Assertion failures will only happen in the debug mode tests;
		// release mode tests won't trip this exception.
		try {
			set_throw_on_next_assertion_failure();
			r1_lkbinfo.atom1_derivs();
			TS_ASSERT( false ); // this line should never execute
		} catch ( utility::excn::Exception & e ) {
		}
#endif

		// Now let's try again, but where we've run setup_for_derivatives, and
		// so the water derivative matrices have been allocated & computed.

		ResSingleMinimizationData mindat;
		core::scoring::methods::EnergyMethodOptions options;
		core::scoring::lkball::LK_BallEnergy lkbE( options );

		basic::datacache::BasicDataCache & r1_cache( pose.residue_data(1) );
		lkbE.setup_for_derivatives_for_residue( res1, pose, sfxn, mindat, r1_cache );

		// residue 1's data cache should now include derivative info
		try {
			set_throw_on_next_assertion_failure();
			utility::vector1< WaterDerivMatrix > at2_derivs = r1_lkbinfo.atom2_derivs();
			TS_ASSERT_EQUALS( at2_derivs.size(), n_waters_actual );
		} catch ( utility::excn::Exception & e ) {
			TS_ASSERT( false ); // this line should never execute
		}

	}


	void test_lk_ball_rpe_matches_bbbb_bbsc_scscE()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( core::scoring::lk_ball, 1.0 );
		core::scoring::methods::EnergyMethodOptions options;

		core::scoring::lkball::LK_BallEnergy lkbE( options );
		TS_ASSERT( lkbE.divides_backbone_and_sidechain_energetics() );
		sfxn( pose );

		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			for ( Size jj = ii+1; jj <= pose.total_residue(); ++jj ) {
				EnergyMap emap_rpe, emap_divided;
				lkbE.residue_pair_energy( pose.residue( ii ), pose.residue( jj ), pose, sfxn, emap_rpe );

				lkbE.backbone_backbone_energy( pose.residue( ii ), pose.residue( jj ), pose, sfxn, emap_divided );
				lkbE.backbone_sidechain_energy( pose.residue( ii ), pose.residue( jj ), pose, sfxn, emap_divided );
				lkbE.backbone_sidechain_energy( pose.residue( jj ), pose.residue( ii ), pose, sfxn, emap_divided );
				lkbE.sidechain_sidechain_energy( pose.residue( ii ), pose.residue( jj ), pose, sfxn, emap_divided );

				for ( Size kk = 1; kk <= lkbE.score_types().size(); ++kk ) {
					ScoreType kkst = lkbE.score_types()[ kk ];
					TS_ASSERT_DELTA( emap_rpe[ kkst ], emap_divided[ kkst ], 1e-6 );
				}
			}
		}
	}

	void test_lk_ball_bridge_rpe_matches_bbbb_bbsc_scscE()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( core::scoring::lk_ball, 1.0 );
		sfxn.set_weight( core::scoring::lk_ball_bridge, 1.0 );
		core::scoring::methods::EnergyMethodOptions options;

		core::scoring::lkball::LK_BallEnergy lkbE( options );
		TS_ASSERT( lkbE.divides_backbone_and_sidechain_energetics() );
		sfxn( pose );

		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			for ( Size jj = ii+1; jj <= pose.total_residue(); ++jj ) {
				EnergyMap emap_rpe, emap_divided;
				lkbE.residue_pair_energy( pose.residue( ii ), pose.residue( jj ), pose, sfxn, emap_rpe );

				lkbE.backbone_backbone_energy( pose.residue( ii ), pose.residue( jj ), pose, sfxn, emap_divided );
				lkbE.backbone_sidechain_energy( pose.residue( ii ), pose.residue( jj ), pose, sfxn, emap_divided );
				lkbE.backbone_sidechain_energy( pose.residue( jj ), pose.residue( ii ), pose, sfxn, emap_divided );
				lkbE.sidechain_sidechain_energy( pose.residue( ii ), pose.residue( jj ), pose, sfxn, emap_divided );

				for ( Size kk = 1; kk <= lkbE.score_types().size(); ++kk ) {
					ScoreType kkst = lkbE.score_types()[ kk ];
					TS_ASSERT_DELTA( emap_rpe[ kkst ], emap_divided[ kkst ], 1e-6 );
				}
			}
		}
	}

	void test_start_func_matches_start_score_w_full_bbflex()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		core::scoring::methods::EnergyMethodOptions options;
		options.analytic_etable_evaluation( false );
		sfxn.set_energy_method_options( options );

		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( lk_ball, 0.625 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );

		/*Size before_precision = std::cout.precision();
		std::cout.precision( 16 );
		core::Real score = sfxn.score( pose );
		std::cout << "MINE: " << score << std::endl;
		std::cout.precision( before_precision );*/

		adv.validate_start_func_matches_start_score( -7.70049544775657, false );
	}

	void test_lkball_numeric_deriv_check()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		core::scoring::methods::EnergyMethodOptions options;
		options.analytic_etable_evaluation( false );
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( lk_ball, 0.625 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv( pose, sfxn, movemap );

		adv.simple_deriv_check( true, 5e-3 );
	}


	/// @brief Create a move map that has the central residue changing, while everything
	/// else is held fixed; then the domain map will have color 1 for residues 1-9, and
	/// color 2 for residues 11-20.
	void test_start_func_matches_start_score_w_partial_bbflex()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		core::scoring::methods::EnergyMethodOptions options;
		options.analytic_etable_evaluation( false );
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( lk_ball, 0.625 );
		kinematics::MoveMap movemap( create_trpcage_movemap_to_allow_bb10_freedom() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );

		/*Size before_precision = std::cout.precision();
		std::cout.precision( 16 );
		core::Real score = sfxn.score( pose );
		std::cout << "MINE: " << score << std::endl;
		std::cout.precision( before_precision );*/

		adv.validate_start_func_matches_start_score( -7.70049544775657, false );

		adv.simple_deriv_check( true, 5e-3 );
	}

	void test_nblist_auto_update_derivatives_w_full_bb_flex()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		core::scoring::methods::EnergyMethodOptions options;
		options.analytic_etable_evaluation( false );
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( lk_ball, 0.625 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.set_nblist_auto_update( true );
		adv.validate_start_func_matches_start_score( -7.70049544775657, false );

		adv.simple_deriv_check( true, 5e-3 );
	}

	void test_nblist_auto_update_derivatives_w_partial_bb_flex()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		core::scoring::methods::EnergyMethodOptions options;
		options.analytic_etable_evaluation( false );
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( lk_ball, 0.625 );
		kinematics::MoveMap movemap( create_trpcage_movemap_to_allow_bb10_freedom() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.set_nblist_auto_update( true );
		adv.validate_start_func_matches_start_score( -7.70049544775657, false );

		adv.simple_deriv_check( true, 5e-3 );
	}

};
