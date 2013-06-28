// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/hbonds/HBondEnergy.cxxtest.hh
/// @brief  Test the score and derivative evaluation for the HBondEnergy class
/// @author Andrew Leaver-Fay


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/UTracer.hh>
#include <test/util/deriv_funcs.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Package headers
// AUTO-REMOVED #include <core/scoring/hbonds/constants.hh>
// AUTO-REMOVED #include <core/scoring/hbonds/hbonds_geom.hh>
// AUTO-REMOVED #include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/HBondEnergy.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
// AUTO-REMOVED #include <core/scoring/hbonds/HBondTypeManager.hh>

#include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
// AUTO-REMOVED #include <core/optimization/MinimizerMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>

// AUTO-REMOVED #include <basic/Tracer.hh>

// AUTO-REMOVED #include <core/id/DOF_ID.hh>
// AUTO-REMOVED #include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionInfo.hh>
#include <core/scoring/EnergyGraph.hh>

#include <core/import_pose/import_pose.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>

// Utility headers
#include <utility/vector1.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/optimization/AtomTreeMinimizer.fwd.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/optimization/NumericalDerivCheckResult.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/SecondaryStructureWeights.hh>
#include <core/scoring/hbonds/HBondOptions.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/mm/MMBondAngleResidueTypeParamSet.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/xyzVector.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

//Auto using namespaces
namespace std { } using namespace std; // AUTO USING NS
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end


using namespace core;
using namespace conformation;
using namespace scoring;
using namespace hbonds;


class HBondEnergyTests : public CxxTest::TestSuite {

public:
  void setUp() {
		core_init();
	}

	void tearDown(){}

	/// @brief Setup for minimization using a move map that says "minimize all bb and sc torsions"
	/// Make sure that start_score matches start_func.
	void test_hbonds_deriv_check_w_full_torsional_flexibility()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		//sfxn.set_weight( fa_atr, 0.5 ); /// just test hbonds; ignore etable
		//sfxn.set_weight( fa_rep, 0.5 );
		sfxn.set_weight( hbond_sr_bb, 0.75 );
		sfxn.set_weight( hbond_lr_bb, 0.5  );
		sfxn.set_weight( hbond_bb_sc, 1.25 );
		sfxn.set_weight( hbond_sc, 1.125 );

		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		AtomDerivValidator adv( pose, sfxn, movemap );

		/// This call runs a numeric deriv check on all the free dofs in the system and makes sure
		/// that the analytic norm matches the numeric norm to 1e-4.
		adv.simple_deriv_check( true, 1e-6 );
	}

	void test_hbonds_w_sp2term_deriv_check_w_full_torsional_flexibility()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		EnergyMethodOptionsOP emopts( new EnergyMethodOptions( sfxn.energy_method_options() ));
		emopts->hbond_options().use_sp2_chi_penalty( true );
		sfxn.set_energy_method_options( *emopts );
		//sfxn.set_weight( fa_atr, 0.5 ); /// just test hbonds; ignore etable
		//sfxn.set_weight( fa_rep, 0.5 );
		sfxn.set_weight( hbond_sr_bb, 0.75 );
		sfxn.set_weight( hbond_lr_bb, 0.5  );
		sfxn.set_weight( hbond_bb_sc, 1.25 );
		sfxn.set_weight( hbond_sc, 1.125 );

		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		AtomDerivValidator adv( pose, sfxn, movemap );

		/// This call runs a numeric deriv check on all the free dofs in the system and makes sure
		/// that the analytic norm matches the numeric norm to 1e-4.
		adv.simple_deriv_check( true, 1e-6 );

	}

	/// @brief This cannot be made a unit test, because it is very unlikely that multiple machines
	/// would minimize to the same structure and score.  There is simply too much numerical noise
	/// in the process of minimization
	void dont_test_atom_tree_minimize_with_hbonds_and_etable()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( hbond_sr_bb, 0.75 );
		sfxn.set_weight( hbond_lr_bb, 0.5  );
		sfxn.set_weight( hbond_bb_sc, 1.25 );
		sfxn.set_weight( hbond_sc, 1.125 );

		AtomTreeMinimizer minimizer;
		Real start_score = sfxn(pose);
		std::cout.precision( 16 );
		//std::cout << "start score: " << sfxn(pose) << std::endl;
		TS_ASSERT_DELTA( -33.67974148760877, start_score, 1e-12 );


		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, true, false );

		minimizer.run( pose, movemap, sfxn, min_options );

		Real end_score = sfxn(pose);
		//std::cout << "end score: " << sfxn(pose) << std::endl;
		TS_ASSERT_DELTA( -36.48274046926093, end_score, 1e-12 );
	}

	/// @brief Create a move map that has the central residue changing, while everything
	/// else is held fixed; then the domain map will have color 1 for residues 1-9, and
	/// color 2 for residues 11-20.
	void test_hbonds_deriv_check_w_partial_torsional_flexibility()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		EnergyMethodOptionsOP emopts( new EnergyMethodOptions( sfxn.energy_method_options() ));
		emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );
		sfxn.set_energy_method_options( *emopts );

		sfxn.set_weight( hbond_sr_bb, 0.75 );
		sfxn.set_weight( hbond_lr_bb, 0.5  );
		sfxn.set_weight( hbond_bb_sc, 1.25 );
		sfxn.set_weight( hbond_sc, 1.125 );

		kinematics::MoveMap movemap;
		movemap.set_bb( 10, true );
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-7 );
	}

	/// @brief Create a move map that has the central residue changing, while everything
	/// else is held fixed; then the domain map will have color 1 for residues 1-9, and
	/// color 2 for residues 11-20.
	void test_hbonds_w_sp2term_deriv_check_w_partial_torsional_flexibility()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		EnergyMethodOptionsOP emopts( new EnergyMethodOptions( sfxn.energy_method_options() ));
		emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );
		emopts->hbond_options().use_sp2_chi_penalty( true );
		sfxn.set_energy_method_options( *emopts );

		sfxn.set_weight( hbond_sr_bb, 0.75 );
		sfxn.set_weight( hbond_lr_bb, 0.5  );
		sfxn.set_weight( hbond_bb_sc, 1.25 );
		sfxn.set_weight( hbond_sc, 1.125 );

		kinematics::MoveMap movemap;
		movemap.set_bb( 10, true );
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-7 );
	}

	void dont_test_atom_tree_minimize_with_hbonds_and_etable2()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.5 );
		sfxn.set_weight( hbond_sr_bb, 0.75 );
		sfxn.set_weight( hbond_lr_bb, 0.5  );
		sfxn.set_weight( hbond_bb_sc, 1.25 );
		sfxn.set_weight( hbond_sc, 1.125 );

		AtomTreeMinimizer minimizer;
		Real start_score = sfxn(pose);
		std::cout.precision( 16 );
		//std::cout << "start score: " << sfxn(pose) << std::endl;
		TS_ASSERT_DELTA( -33.67974148760877, start_score, 1e-12 );


		kinematics::MoveMap movemap;
		movemap.set_bb( 10, true );

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, false, false );

		minimizer.run( pose, movemap, sfxn, min_options );

		Real end_score = sfxn(pose);
		//std::cout << "end score: " << sfxn(pose) << std::endl;
		TS_ASSERT_DELTA( -33.88566264995315, end_score, 1e-12 );
	}

	void initialize_lj_hbond_sfxn( core::scoring::ScoreFunction & sfxn )
	{
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.5 );
		sfxn.set_weight( hbond_sr_bb, 0.75 );
		sfxn.set_weight( hbond_lr_bb, 0.5  );
		sfxn.set_weight( hbond_bb_sc, 1.25 );
		sfxn.set_weight( hbond_sc, 1.125 );
	}

	void create_excludable_bb_sc_hbond( core::pose::Pose & pose )
	{
		using namespace core::id;

		core::import_pose::pose_from_pdb( pose, "core/scoring/hbonds/1ubq_23_to_34.pdb" ); // read in the helix from ubiquitin
		core::scoring::ScoreFunction sfxn;
		initialize_lj_hbond_sfxn( sfxn );

		std::string resfile_mutate_l26s( "NATRO\nstart\n30 _ PIKAA S\n" );
		core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( pose );
		core::pack::task::parse_resfile_string( pose, *task, resfile_mutate_l26s );
		core::pack::pack_rotamers( pose, sfxn, task );

		TS_ASSERT( pose.residue_type( 8 ).aa() == chemical::aa_ser );

		core::Size res8_HG_ind = pose.residue_type( 8 ).atom_index( "HG" );
		core::Size res4_O_ind = pose.residue_type( 4 ).atom_index( "O" );

		//Real closest_chi1( 0 ), closest_chi2( 0 ), mind2( -1 );
		//for ( core::Real chi1= -80; chi1 < -40; chi1 += 0.5 ) {
		//	for ( core::Real chi2 = 60; chi2 < 80; chi2 += 0.5 ) {
		//		pose.set_chi( 1, 8, chi1 );
		//		pose.set_chi( 2, 8, chi2 );
		//		core::Real d2 = pose.xyz( AtomID( res8_HG_ind, 8 )).distance_squared( pose.xyz( AtomID( res4_O_ind, 4 )) );
		//		if ( mind2 < 0 || d2 < mind2 ) {
		//			closest_chi1 = chi1;
		//			closest_chi2 = chi2;
		//			mind2 = d2;
		//		}
		//	}
		//}
		//std::cout << "closest res8_HG res4_0 distance " << std::sqrt( mind2 ) << " with chi1= " << closest_chi1 << " and chi2= " << closest_chi2 << std::endl;


		// set chi for residue 4 to form hbond
		pose.set_chi( 1, 8, -40.5 );
		pose.set_chi( 2, 8, 60 );

		// HG on residue 30 should (residue 8) should be ~2.5 A from the backbone O on residue 26 (residue 4 )
		TS_ASSERT( pose.xyz( AtomID( res8_HG_ind, 8 )).distance( pose.xyz( AtomID( res4_O_ind, 4 )) ) < 2.6 );

		pose.dump_pdb( "excludable_bb_sc_hbond.pdb" );
	}

	void test_hbonds_w_excluded_bb_sc_interactions() {
		using namespace core::pose;
		using namespace core::scoring;

		Pose pose;
		create_excludable_bb_sc_hbond( pose );
		ScoreFunction sfxn;
		initialize_lj_hbond_sfxn( sfxn );
		sfxn( pose );

		{ // scope
		core::graph::Edge const * e4_8 = pose.energies().energy_graph().find_edge( 4, 8 );
		TS_ASSERT( e4_8 );
		if ( ! e4_8 ) return;

		EnergyEdge const * ee4_8 = dynamic_cast< EnergyEdge const * > (e4_8);
		TS_ASSERT( ee4_8 );
		if ( ! ee4_8 ) return;

		// there should be no hydrogen bond detected, since residue 4 is already participating in a bb/bb hbond.
		TS_ASSERT( (*ee4_8)[ hbond_bb_sc ] == 0.0 );
		}

		ScoreFunction sfxn2;
		methods::EnergyMethodOptions enmethopts = sfxn2.energy_method_options();
		enmethopts.hbond_options().bb_donor_acceptor_check( false ); // <-- disable the bb/sc exclusion rule.
		sfxn2.set_energy_method_options( enmethopts );
		initialize_lj_hbond_sfxn( sfxn2 );
		TS_ASSERT( pose.energies().get_scorefxn_info() != (* sfxn2.info()) );

		sfxn2( pose );

		{ // scope
		core::graph::Edge const * e4_8 = pose.energies().energy_graph().find_edge( 4, 8 );
		TS_ASSERT( e4_8 );
		if ( ! e4_8 ) return;

		EnergyEdge const * ee4_8 = dynamic_cast< EnergyEdge const * > (e4_8);
		TS_ASSERT( ee4_8 );
		if ( ! ee4_8 ) return;

		// turning off the bb/sc exclusion rule should produce an hbond.
		TS_ASSERT( (*ee4_8)[ hbond_bb_sc ] != 0.0 );
		}

		HBondEnergy hbe_w_bbsc_exclusion(   sfxn.energy_method_options().hbond_options() );
		HBondEnergy hbe_wo_bbsc_exclusion( sfxn2.energy_method_options().hbond_options() );

		// Let's check the backbone_sidechain_energy() function to make sure the exclusion rule is properly applied.
		EnergyMap emap;
		sfxn( pose );
		hbe_w_bbsc_exclusion.backbone_sidechain_energy( pose.residue(4), pose.residue(8), pose, sfxn, emap );
		TS_ASSERT( emap[ hbond_bb_sc ] == 0.0 );

		emap.zero();
		sfxn2( pose );
		hbe_wo_bbsc_exclusion.backbone_sidechain_energy( pose.residue(4), pose.residue(8), pose, sfxn, emap );
		TS_ASSERT( emap[ hbond_bb_sc ] != 0.0 );

	}

	void verify_hbond_trie_vs_trie_calculation( ScoreFunction const & sfxn )
	{
		using namespace core::pose;
		using namespace core::scoring;

		Pose pose;
		create_excludable_bb_sc_hbond( pose );
		sfxn( pose );

		core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( pose );
		// design residues 4 and 8. Turn everything else off.
		for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			if ( ii == 4 || ii == 8 ) continue;
			task->nonconst_residue_task(ii).prevent_repacking();
		}

		task->nonconst_residue_task(4).and_extrachi_cutoff( 1 );
		task->nonconst_residue_task(8).and_extrachi_cutoff( 1 );

		//pack_scorefxn_pose_handshake( pose, sfxn);

		//replace this with RotSetsFactory
		core::pack::rotamer_set::RotamerSetsOP rotsets( new core::pack::rotamer_set::RotamerSets() );
		core::pack::interaction_graph::InteractionGraphBaseOP ig = NULL;

		// create the rotamer sets so that we can then compute the trie-vs-trie energies
		core::pack::pack_rotamers_setup( pose, sfxn, task, rotsets, ig );

		core::Size res4_nrots = rotsets->rotamer_set_for_residue( 4 )->num_rotamers();
		core::Size res8_nrots = rotsets->rotamer_set_for_residue( 8 )->num_rotamers();

		ObjexxFCL::FArray2D< core::PackerEnergy > rpe_table( res8_nrots, res4_nrots, core::PackerEnergy( 0.0 ) );
		HBondEnergy hbe( sfxn.energy_method_options().hbond_options() );
		hbe.evaluate_rotamer_pair_energies(
			*rotsets->rotamer_set_for_residue( 4 ),
			*rotsets->rotamer_set_for_residue( 8 ),
			pose, sfxn, sfxn.weights(), rpe_table );

		EnergyMap emap;
		for ( Size ii = 1; ii <= res4_nrots; ++ii ) {
			for ( Size jj = 1; jj <= res8_nrots; ++jj ) {
				emap.zero();
				hbe.residue_pair_energy(
					*rotsets->rotamer_set_for_residue(4)->rotamer( ii ),
					*rotsets->rotamer_set_for_residue(8)->rotamer( jj ),
					pose, sfxn, emap );
				core::PackerEnergy score = emap.dot( sfxn.weights() );
				TS_ASSERT_DELTA( score, rpe_table( jj, ii ), 1e-6 );
			}
		}


	}

	/// Ensure that the rotamer pair energy evaluation works correctly when using the trie-vs-trie algorithm
	void test_hbonds_trie_vs_trie_standard() {

		ScoreFunction sfxn;
		methods::EnergyMethodOptions enmethopts = sfxn.energy_method_options();
		enmethopts.hbond_options().decompose_bb_hb_into_pair_energies( true ); // calculate bb/bb hbonds in residue_pair_energy calls.
		sfxn.set_energy_method_options( enmethopts );

		initialize_lj_hbond_sfxn( sfxn );

		verify_hbond_trie_vs_trie_calculation( sfxn );
	}

	/// Ensure that the rotamer pair energy evaluation works correctly when using the trie-vs-trie algorithm
	void test_hbonds_trie_vs_trie_no_bbsc_exclusion() {

		ScoreFunction sfxn;
		methods::EnergyMethodOptions enmethopts = sfxn.energy_method_options();
		enmethopts.hbond_options().decompose_bb_hb_into_pair_energies( true ); // calculate bb/bb hbonds in residue_pair_energy calls.
		enmethopts.hbond_options().bb_donor_acceptor_check( false ); // <-- disable the bb/sc exclusion rule.

		sfxn.set_energy_method_options( enmethopts );

		initialize_lj_hbond_sfxn( sfxn );

		verify_hbond_trie_vs_trie_calculation( sfxn );
	}

	/// Ensure that the rotamer pair energy evaluation works correctly when using the trie-vs-trie algorithm
	void test_hbonds_trie_vs_trie_sp2() {

		ScoreFunction sfxn;
		methods::EnergyMethodOptions enmethopts = sfxn.energy_method_options();
		enmethopts.hbond_options().decompose_bb_hb_into_pair_energies( true ); // calculate bb/bb hbonds in residue_pair_energy calls.
		enmethopts.hbond_options().use_sp2_chi_penalty( true ); // activate the sp2 potential

		sfxn.set_energy_method_options( enmethopts );

		initialize_lj_hbond_sfxn( sfxn );

		verify_hbond_trie_vs_trie_calculation( sfxn );
	}

};
