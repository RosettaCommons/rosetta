// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/etable/EtableEnergy.cxxtest.hh
/// @brief  Unit tests for the lennard-jones and EEF1 solvation model.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/deriv_funcs.hh>

#include <core/types.hh>

// Unit headers
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/etrie/EtableAtom.hh>
#include <core/scoring/etable/etrie/CountPairData_1_2.hh>
#include <core/scoring/etable/EtableOptions.hh>

// Package headers
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/trie/RotamerTrie.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

#include <core/graph/Graph.hh>
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

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION



using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.etable.EtableEnergy.cxxtest");

// using declarations
using namespace std;
using namespace core;
using namespace scoring;
using namespace etable;


///////////////////////////////////////////////////////////////////////////
/// @name EtableEnergyTest
/// @brief: Test the functionality of the EtableEnergy class
///////////////////////////////////////////////////////////////////////////
class EtableEnergyTests : public CxxTest::TestSuite {

public:

	void setUp() {
		/// The analytic version operates just fine in the presence of the table-based etable, but
		/// the table-based version won't work unless the analytic_etable_evaluation flag is false
		/// at the time the FA_STANDARD Etable is constructed.
		core_init_with_additional_options( "-analytic_etable_evaluation false" );
	}

	void tearDown() {}

	void test_default_enmethopts_signal_analytic_etable() {
		core_init(); // reinitialize options system without the "-analytic_etable_evaluation false" override in setUp above
		core::scoring::methods::EnergyMethodOptions opts;
		TS_ASSERT( opts.analytic_etable_evaluation() );
	}

	void test_eval_residue_pair_energy_table_version()
	{
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;

		Pose pose = create_trpcage_ideal_pose();
		EnergyMethodOptions options;
		options.analytic_etable_evaluation( false );

		TableLookupEtableEnergy etab_energy( *( ScoringManager::get_instance()->etable( options.etable_type() ).lock()), options );

		EnergyMap emap;
		ScoreFunction sfxn;
		etab_energy.residue_pair_energy( pose.residue( 4 ), pose.residue( 5 ), pose, sfxn, emap );

		//int const old_precision( TR.precision() );
		//TR.precision( 16 );
		//TR << "etable energy, trpcage 4 and 5: " << emap[ fa_atr ] << " " << emap[ fa_rep ] << " " << emap[ fa_sol ] << std::endl;
		//TR.precision( old_precision );

		//core.scoring.etable.EtableEnergy.cxxtest: etable energy, trpcage 4 and 5: -0.6079854186280232 0.0778507865338765 0.4275630856840362
		TS_ASSERT_DELTA( emap[ fa_atr ], -0.6079854186280232, 1e-12 );
		TS_ASSERT_DELTA( emap[ fa_rep ],  0.0778507865338765, 1e-12 );
		TS_ASSERT_DELTA( emap[ fa_sol ],  0.4041439359347673, 1e-12 );
	}

	void test_eval_residue_pair_energy_analytic_version()
	{
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;

		Pose pose = create_trpcage_ideal_pose();
		EnergyMethodOptions options;
		options.analytic_etable_evaluation( true );

		AnalyticEtableEnergy etab_energy( *( ScoringManager::get_instance()->etable( options.etable_type() ).lock()), options );

		EnergyMap emap;
		ScoreFunction sfxn;
		etab_energy.residue_pair_energy( pose.residue( 4 ), pose.residue( 5 ), pose, sfxn, emap );

		//int const old_precision( TR.precision() );
		//TR.precision( 16 );
		//TR << "etable energy, trpcage 4 and 5: " << emap[ fa_atr ] << " " << emap[ fa_rep ] << " " << emap[ fa_sol ] << std::endl;
		//TR.precision( old_precision );

		//core.scoring.etable.EtableEnergy.cxxtest: etable energy, trpcage 4 and 5: -0.6079854186280232 0.0778507865338765 0.4275630856840362
		TS_ASSERT_DELTA( emap[ fa_atr ], -0.607982997381734, 1e-12 );
		TS_ASSERT_DELTA( emap[ fa_rep ],  0.0778346123883227, 1e-12 );
		TS_ASSERT_DELTA( emap[ fa_sol ],  0.40414285103945, 1e-12 );
	}

	void test_eval_residue_pair_energy_w_minimization_data()
	{
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;

		Pose pose = create_trpcage_ideal_pose();
		EnergyMethodOptions options; // default is fine

		TableLookupEtableEnergy etab_energy( *( ScoringManager::get_instance()->etable( options.etable_type() ).lock()), options );

		EnergyMap emap;
		ScoreFunction sfxn;
		etab_energy.residue_pair_energy( pose.residue( 4 ), pose.residue( 5 ), pose, sfxn, emap );

		optimization::MinimizerMap minmap;
		EnergyMap emap2;
		ResSingleMinimizationData r4dat, r5dat;
		ResPairMinimizationData min_data;
		etab_energy.setup_for_minimizing_for_residue( pose.residue( 4 ), pose, sfxn, minmap, r4dat );
		etab_energy.setup_for_minimizing_for_residue( pose.residue( 5 ), pose, sfxn, minmap, r5dat );
		etab_energy.setup_for_minimizing_for_residue_pair( pose.residue( 4 ), pose.residue( 5 ), pose, sfxn,  minmap, r4dat, r5dat, min_data );
		etab_energy.residue_pair_energy_ext( pose.residue( 4 ), pose.residue( 5 ), min_data, pose, sfxn, emap2 );

		TS_ASSERT_DELTA( emap[ fa_atr ], emap2[ fa_atr ], 1e-12 );
		TS_ASSERT_DELTA( emap[ fa_rep ], emap2[ fa_rep ], 1e-12 );
		TS_ASSERT_DELTA( emap[ fa_sol ], emap2[ fa_sol ], 1e-12 );
	}

	void test_eval_intra_residue_energy_w_minimization_data_table_version()
	{
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;

		Pose pose = create_trpcage_ideal_pose();
		EnergyMethodOptions options; // default is fine

		TableLookupEtableEnergy etab_energy( *( ScoringManager::get_instance()->etable( options.etable_type() ).lock()), options, true );

		EnergyMap emap;
		ScoreFunction sfxn;
		sfxn.set_weight( fa_intra_atr, 0.1 ); // enable intra-residue energy evaluation; skipped if all intra-residue weights are 0
		etab_energy.eval_intrares_energy( pose.residue( 4 ), pose, sfxn, emap );

		optimization::MinimizerMap minmap;
		EnergyMap emap2;
		ResSingleMinimizationData min_data;
		etab_energy.setup_for_minimizing_for_residue( pose.residue( 4 ), pose, sfxn,  minmap, min_data );
		etab_energy.eval_intrares_energy_ext( pose.residue( 4 ), min_data, pose, sfxn, emap2 );

		TS_ASSERT_DELTA( emap[ fa_intra_atr ], emap2[ fa_intra_atr ], 1e-12 );
		TS_ASSERT_DELTA( emap[ fa_intra_rep ], emap2[ fa_intra_rep ], 1e-12 );
		TS_ASSERT_DELTA( emap[ fa_intra_sol ], emap2[ fa_intra_sol ], 1e-12 );

	}

	void test_eval_intra_residue_energy_w_minimization_data_analytic_version()
	{
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;

		Pose pose = create_trpcage_ideal_pose();
		EnergyMethodOptions options; // default is fine

		AnalyticEtableEnergy etab_energy( *( ScoringManager::get_instance()->etable( options.etable_type() ).lock()), options, true );

		EnergyMap emap;
		ScoreFunction sfxn;
		sfxn.set_weight( fa_intra_atr, 0.1 ); // enable intra-residue energy evaluation; skipped if all intra-residue weights are 0
		etab_energy.eval_intrares_energy( pose.residue( 4 ), pose, sfxn, emap );

		optimization::MinimizerMap minmap;
		EnergyMap emap2;
		ResSingleMinimizationData min_data;
		etab_energy.setup_for_minimizing_for_residue( pose.residue( 4 ), pose, sfxn,  minmap, min_data );
		etab_energy.eval_intrares_energy_ext( pose.residue( 4 ), min_data, pose, sfxn, emap2 );

		TS_ASSERT_DELTA( emap[ fa_intra_atr ], emap2[ fa_intra_atr ], 1e-12 );
		TS_ASSERT_DELTA( emap[ fa_intra_rep ], emap2[ fa_intra_rep ], 1e-12 );
		TS_ASSERT_DELTA( emap[ fa_intra_sol ], emap2[ fa_intra_sol ], 1e-12 );

	}

	/*void dont_test_atom_deriv_old_vs_new()
	{
	using namespace core;
	using namespace core::graph;
	using namespace core::pose;
	using namespace core::scoring;
	using namespace core::scoring::etable;
	using namespace core::scoring::methods;

	Pose pose = create_trpcage_ideal_pose();
	EnergyMethodOptions options; // default is fine

	EtableEnergy etab_energy( *( ScoringManager::get_instance()->etable( options.etable_type() )), options );

	EnergyMap emap;
	ScoreFunction sfxn;
	sfxn.set_weight( fa_atr, 0.5 );
	sfxn.set_weight( fa_rep, 0.25 );
	sfxn.set_weight( fa_sol, 0.125 );

	//sfxn.set_weight( fa_intra_rep, 0.02 );

	sfxn( pose );
	optimization::MinimizerMap minmap;
	kinematics::MoveMap movemap;
	movemap.set_bb( true );
	movemap.set_chi( true );
	minmap.setup( pose, movemap );

	MinimizationGraph g( pose.size() );
	g.copy_connectivity( pose.energies().energy_graph() );
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
	etab_energy.setup_for_minimizing_for_residue(
	pose.residue( ii ), pose, sfxn, minmap, g.get_minimization_node( ii )->res_min_data() );
	}
	for ( Graph::EdgeListIter edge_iter = g.edge_list_begin(),
	edge_iter_end = g.edge_list_end(); edge_iter != edge_iter_end;
	++edge_iter ) {
	Size const node1 = (*edge_iter)->get_first_node_ind();
	Size const node2 = (*edge_iter)->get_second_node_ind();
	MinimizationEdge * minedge = dynamic_cast< MinimizationEdge * > ( *edge_iter );
	//std::cout << "Setting up minimization data for pair " << node1 << " " << node2 << std::endl;
	etab_energy.setup_for_minimizing_for_residue_pair(
	pose.residue( node1 ), pose.residue( node2 ), pose, sfxn,  minmap,
	g.get_minimization_node( node1 )->res_min_data(),
	g.get_minimization_node( node2 )->res_min_data(),
	minedge->res_pair_min_data() );
	}
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
	MinimizationNode * iinode = g.get_minimization_node( ii );
	etab_energy.setup_for_minimizing_for_residue( pose.residue( ii ), pose, sfxn, minmap, iinode->res_min_data() );
	}

	/// Let's measure the derivatie for the CD2 atom on residue 6
	id::AtomID cd2( pose.residue(6).atom_index("CD2"), 6 );

	pose.energies().set_use_nblist( pose, minmap.domain_map(), false ); // false = autoupdate
	sfxn.setup_for_minimizing( pose, minmap );
	sfxn.setup_for_derivatives( pose );

	Vector goldF1(0.0), goldF2(0.0);
	sfxn.eval_atom_derivative( cd2, pose, minmap.domain_map(), goldF1, goldF2 );


	Vector decompF1(0.0), decompF2(0.0);
	EnergyMap weights;
	weights[ fa_atr ] = 0.5;
	weights[ fa_rep ] = 0.25;
	weights[ fa_sol ] = 0.125;
	//weights[ fa_intra_rep ] = 0.02;
	for ( Graph::EdgeListConstIter edge_iter = g.get_node( 6 )->const_edge_list_begin(),
	edge_iter_end = g.get_node( 6 )->const_edge_list_end(); edge_iter != edge_iter_end; ++edge_iter ) {
	MinimizationEdge const * minedge = dynamic_cast< MinimizationEdge const * > ( *edge_iter );
	Size const node2 = (*edge_iter)->get_other_ind( 6 );
	//std::cout << "Evaluating derivative for pair " << (*edge_iter)->get_first_node_ind() << " " << node2 << std::endl;
	etab_energy.eval_atom_derivative_for_residue_pair(
	cd2.atomno(), pose.residue( 6 ), pose.residue( node2 ),
	g.get_minimization_node( 6 )->res_min_data(), g.get_minimization_node( node2 )->res_min_data(),
	minedge->res_pair_min_data(), pose, minmap.domain_map(),
	sfxn, weights, decompF1, decompF2 );
	}
	etab_energy.eval_intrares_atom_derivative( cd2.atomno(), pose.residue( 6 ),
	g.get_minimization_node(6)->res_min_data(), pose, minmap.domain_map(),
	sfxn, weights, decompF1, decompF2 );

	TS_ASSERT_DELTA( goldF1.x(), decompF1.x(), 1e-12 );
	TS_ASSERT_DELTA( goldF1.y(), decompF1.y(), 1e-12 );
	TS_ASSERT_DELTA( goldF1.z(), decompF1.z(), 1e-12 );
	TS_ASSERT_DELTA( goldF2.x(), decompF2.x(), 1e-12 );
	TS_ASSERT_DELTA( goldF2.y(), decompF2.y(), 1e-12 );
	TS_ASSERT_DELTA( goldF2.z(), decompF2.z(), 1e-12 );
	}*/

	void test_start_func_matches_start_score_w_full_bbflex_table_version()
	{

		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		core::scoring::methods::EnergyMethodOptions options;
		options.analytic_etable_evaluation( false );
		sfxn.set_energy_method_options( options );

		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.validate_start_func_matches_start_score( -26.03949788109572 );
	}

	void test_start_func_matches_start_score_w_full_bbflex_analytic_version()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		core::scoring::methods::EnergyMethodOptions options;
		options.analytic_etable_evaluation( false );
		sfxn.set_energy_method_options( options );

		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.validate_start_func_matches_start_score( -26.03949788109572 );
	}

	void test_etab_numeric_deriv_check_table_version()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		core::scoring::methods::EnergyMethodOptions options;
		options.analytic_etable_evaluation( false );
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv( pose, sfxn, movemap );

		adv.simple_deriv_check( true, 5e-3 );
	}

	void test_etab_numeric_deriv_check_analytic_version()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		core::scoring::methods::EnergyMethodOptions options;
		options.analytic_etable_evaluation( true );
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv( pose, sfxn, movemap );

		adv.simple_deriv_check( true, 1e-6 );
	}

	/// @brief Create a move map that has the central residue changing, while everything
	/// else is held fixed; then the domain map will have color 1 for residues 1-9, and
	/// color 2 for residues 11-20.
	void test_start_func_matches_start_score_w_partial_bbflex_table_version()
	{

		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		core::scoring::methods::EnergyMethodOptions options;
		options.analytic_etable_evaluation( false );
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		kinematics::MoveMap movemap( create_trpcage_movemap_to_allow_bb10_freedom() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.validate_start_func_matches_start_score( -26.03949788109572 );
	}


	/// @brief Create a move map that has the central residue changing, while everything
	/// else is held fixed; then the domain map will have color 1 for residues 1-9, and
	/// color 2 for residues 11-20.
	void test_start_func_matches_start_score_w_partial_bbflex_analytic_version()
	{

		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		core::scoring::methods::EnergyMethodOptions options;
		options.analytic_etable_evaluation( true );
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		kinematics::MoveMap movemap( create_trpcage_movemap_to_allow_bb10_freedom() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.validate_start_func_matches_start_score( -26.04170800071507 );
	}


	void dont_test_atom_tree_minimize_with_etable_energy2()
	{
		using namespace core;
		using namespace core::graph;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );

		AtomTreeMinimizer minimizer;
		//std::cout << "start score: " << sfxn(pose) << std::endl;
		Real start_score = sfxn(pose);
		TS_ASSERT_DELTA( -25.94299600233435, start_score, 1e-12 );


		kinematics::MoveMap movemap;
		movemap.set_bb( 10, true );

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, false, false );

		minimizer.run( pose, movemap, sfxn, min_options );

		Real end_score = sfxn(pose);
		//std::cout << "end score: " << sfxn(pose) << std::endl;
		TS_ASSERT_DELTA( -26.15005763826676, end_score, 1e-12 );
	}

	void dont_test_atom_tree_minimize_with_intrares_etable_energy()
	{
		using namespace core;
		using namespace core::graph;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );

		AtomTreeMinimizer minimizer;
		//std::cout.precision( 16 );
		//std::cout << "start score: " << sfxn(pose) << std::endl;
		Real start_score = sfxn(pose);
		TS_ASSERT_DELTA( -22.28334391305056, start_score, 1e-12 );


		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, false, false );

		minimizer.run( pose, movemap, sfxn, min_options );

		Real end_score = sfxn(pose);
		//std::cout << "end score: " << sfxn(pose) << std::endl;
		TS_ASSERT_DELTA( -22.82696911663496, end_score, 1e-12 );
	}

	void dont_test_atom_tree_minimize_with_intrares_etable_energy2()
	{
		using namespace core;
		using namespace core::graph;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );

		AtomTreeMinimizer minimizer;
		//std::cout << "start score: " << sfxn(pose) << std::endl;
		Real start_score = sfxn(pose);
		TS_ASSERT_DELTA( -22.28334391305056, start_score, 1e-12 );


		kinematics::MoveMap movemap;
		movemap.set_bb( 10, true );

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, false, false );

		minimizer.run( pose, movemap, sfxn, min_options );

		Real end_score = sfxn(pose);
		//std::cout << "end score: " << sfxn(pose) << std::endl;
		TS_ASSERT_DELTA( -22.49034811587382, end_score, 1e-12 );
	}

	/// @brief Create a move map that has the central residue changing, while everything
	/// else is held fixed; then the domain map will have color 1 for residues 1-9, and
	/// color 2 for residues 11-20.
	void test_start_func_matches_start_score_w_full_bbflex_and_intraresidue_table_version()
	{

		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		core::scoring::methods::EnergyMethodOptions options;
		options.analytic_etable_evaluation( false );
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.validate_start_func_matches_start_score( -22.39185655926555 );
	}

	/// @brief Create a move map that has the central residue changing, while everything
	/// else is held fixed; then the domain map will have color 1 for residues 1-9, and
	/// color 2 for residues 11-20.
	void test_start_func_matches_start_score_w_full_bbflex_and_intraresidue_analytic_version()
	{

		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		core::scoring::methods::EnergyMethodOptions options;
		options.analytic_etable_evaluation( true );
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.validate_start_func_matches_start_score( -22.39669996258689 );
	}

	/// @brief Make sure that the domain map logic inside the ScoreFunction
	/// operates correctly for the intraresidue portions of two-body energies.
	void test_start_func_matches_start_score_w_partial_bbflex_and_intraresidue_table_version()
	{

		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		core::scoring::methods::EnergyMethodOptions options;
		options.analytic_etable_evaluation( false );
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );
		kinematics::MoveMap movemap( create_trpcage_movemap_to_allow_bb10_freedom() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.validate_start_func_matches_start_score( -22.39185655926555 );
	}

	/// @brief Make sure that the domain map logic inside the ScoreFunction
	/// operates correctly for the intraresidue portions of two-body energies.
	void test_start_func_matches_start_score_w_partial_bbflex_and_intraresidue_analytic_version()
	{

		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		core::scoring::methods::EnergyMethodOptions options;
		options.analytic_etable_evaluation( true );
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );
		kinematics::MoveMap movemap( create_trpcage_movemap_to_allow_bb10_freedom() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.validate_start_func_matches_start_score( -22.39669996258689 );
	}


	void test_setup_for_minimizing_with_autoupdate_w_full_bb_flexibility_table_version()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		core::scoring::methods::EnergyMethodOptions options;
		options.analytic_etable_evaluation( false );
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.set_nblist_auto_update( true );
		adv.validate_start_func_matches_start_score( -22.39185655926555 );
	}

	void test_setup_for_minimizing_with_autoupdate_w_full_bb_flexibility_analytic_version()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		core::scoring::methods::EnergyMethodOptions options;
		options.analytic_etable_evaluation( true );
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.set_nblist_auto_update( true );
		adv.validate_start_func_matches_start_score( -22.39669996258689 );
	}

	void test_setup_for_minimizing_with_autoupdate_w_partial_bb_flexibility_table_version()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		core::scoring::methods::EnergyMethodOptions options;
		options.analytic_etable_evaluation( false );
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );
		kinematics::MoveMap movemap( create_trpcage_movemap_to_allow_bb10_freedom() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.set_nblist_auto_update( true );
		adv.validate_start_func_matches_start_score( -22.39185655926555 );

		/*using namespace core;
		using namespace core::graph;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );

		AtomTreeMinimizer minimizer;
		//std::cout.precision( 16 );
		//std::cout << "start score: " << sfxn(pose) << std::endl;
		Real start_score = sfxn(pose);
		TS_ASSERT_DELTA( -22.28334391305056, start_score, 1e-12 );


		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, false, false );
		min_options.nblist_auto_update( true );

		/// BEGIN AtomTreeMinimizer setup block
		MinimizerMap min_map;
		min_map.setup( pose, movemap );

		pose.energies().set_use_nblist( pose, min_map.domain_map(), min_options.nblist_auto_update() );
		sfxn.setup_for_minimizing( pose, min_map );
		/// END AtomTreeMinimizer setup block

		Real start_func = sfxn(pose);
		//std::cout << "start_func: " << start_func << std::endl;
		TS_ASSERT_DELTA( -22.28334391305056, start_func, 1e-12 );*/
	}

	void test_setup_for_minimizing_with_autoupdate_w_partial_bb_flexibility_analytic_version()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		core::scoring::methods::EnergyMethodOptions options;
		options.analytic_etable_evaluation( true );
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );
		kinematics::MoveMap movemap( create_trpcage_movemap_to_allow_bb10_freedom() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.set_nblist_auto_update( true );
		adv.validate_start_func_matches_start_score( -22.39669996258689 );
	}

	void test_etable_derivatives_w_autoupdate_intraresidue_terms_and_full_bb_flex_table_version()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		core::scoring::methods::EnergyMethodOptions options;
		options.analytic_etable_evaluation( false );
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.set_nblist_auto_update( true );

		adv.add_res_for_deriv( 9 );
		adv.add_res_for_deriv( 10 );
		adv.add_res_for_deriv( 11 );

		adv.simple_deriv_check( true, 5e-3 );

	}

	void test_etable_derivatives_w_autoupdate_intraresidue_terms_and_full_bb_flex_analytic_version()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		core::scoring::methods::EnergyMethodOptions options;
		options.analytic_etable_evaluation( true );
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.set_nblist_auto_update( true );

		adv.add_res_for_deriv( 9 );
		adv.add_res_for_deriv( 10 );
		adv.add_res_for_deriv( 11 );

		adv.simple_deriv_check( true, 1e-6 );

	}

	void test_etable_derivatives_w_autoupdate_intraresidue_terms_and_partial_bb_flex_table_version()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		core::scoring::methods::EnergyMethodOptions options;
		options.analytic_etable_evaluation( false );
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );
		kinematics::MoveMap movemap( create_trpcage_movemap_to_allow_bb10_freedom() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.set_nblist_auto_update( true );

		adv.simple_deriv_check( true, 5e-3 );
	}


	void test_etable_derivatives_w_autoupdate_intraresidue_terms_and_partial_bb_flex_analytic_version()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		core::scoring::methods::EnergyMethodOptions options;
		options.analytic_etable_evaluation( true );
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );
		kinematics::MoveMap movemap( create_trpcage_movemap_to_allow_bb10_freedom() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.set_nblist_auto_update( true );

		adv.simple_deriv_check( true, 1e-6 );
	}

	void dont_test_atom_tree_minimize_with_autoupdate()
	{
		using namespace core;
		using namespace core::graph;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );

		AtomTreeMinimizer minimizer;
		//std::cout.precision( 16 );
		//std::cout << "start score: " << sfxn(pose) << std::endl;
		Real start_score = sfxn(pose);
		TS_ASSERT_DELTA( -22.28334391305056, start_score, 1e-12 );


		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, false, false );
		min_options.nblist_auto_update( true );

		minimizer.run( pose, movemap, sfxn, min_options );

		Real end_score = sfxn(pose);
		//std::cout << "end score: " << sfxn(pose) << std::endl;
		TS_ASSERT_DELTA( -22.82696911663496, end_score, 1e-12 );
	}

	void test_setup_for_minimizing_with_autoupdate2_table_version()
	{
		using namespace core;
		using namespace core::graph;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		core::scoring::methods::EnergyMethodOptions options;
		options.analytic_etable_evaluation( false );
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );

		AtomTreeMinimizer minimizer;
		//std::cout.precision( 16 );
		//std::cout << "start score: " << sfxn(pose) << std::endl;
		Real start_score = sfxn(pose);
		TS_ASSERT_DELTA( -22.39185655926555, start_score, 1e-12 );


		kinematics::MoveMap movemap;
		movemap.set_bb( 10, true );

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, false, false );
		min_options.nblist_auto_update( true );

		/// BEGIN AtomTreeMinimizer setup block
		MinimizerMap min_map;
		min_map.setup( pose, movemap );

		pose.energies().set_use_nblist( pose, min_map.domain_map(), min_options.nblist_auto_update() );
		sfxn.setup_for_minimizing( pose, min_map );
		/// END AtomTreeMinimizer setup block

		Real start_func = sfxn(pose);
		//std::cout << "start_func: " << start_func << std::endl;
		TS_ASSERT_DELTA( -22.39185655926555, start_func, 1e-12 );
	}

	void test_setup_for_minimizing_with_autoupdate2_analytic_version()
	{
		using namespace core;
		using namespace core::graph;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		core::scoring::methods::EnergyMethodOptions options;
		options.analytic_etable_evaluation( true );
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );

		AtomTreeMinimizer minimizer;
		//std::cout.precision( 16 );
		//std::cout << "start score: " << sfxn(pose) << std::endl;
		Real start_score = sfxn(pose);
		TS_ASSERT_DELTA( -22.39669996258689, start_score, 1e-12 );


		kinematics::MoveMap movemap;
		movemap.set_bb( 10, true );

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, false, false );
		min_options.nblist_auto_update( true );

		/// BEGIN AtomTreeMinimizer setup block
		MinimizerMap min_map;
		min_map.setup( pose, movemap );

		pose.energies().set_use_nblist( pose, min_map.domain_map(), min_options.nblist_auto_update() );
		sfxn.setup_for_minimizing( pose, min_map );
		/// END AtomTreeMinimizer setup block

		Real start_func = sfxn(pose);
		//std::cout << "start_func: " << start_func << std::endl;
		TS_ASSERT_DELTA( -22.39669996258687, start_func, 1e-12 );
	}

	void dont_test_atom_tree_minimize_with_autoupdate2()
	{
		using namespace core;
		using namespace core::graph;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );

		AtomTreeMinimizer minimizer;
		//std::cout << "start score: " << sfxn(pose) << std::endl;
		Real start_score = sfxn(pose);
		TS_ASSERT_DELTA( -22.28334391305056, start_score, 1e-12 );


		kinematics::MoveMap movemap;
		movemap.set_bb( 10, true );

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, false, false );
		min_options.nblist_auto_update( true );

		minimizer.run( pose, movemap, sfxn, min_options );

		Real end_score = sfxn(pose);
		//std::cout << "end score: " << sfxn(pose) << std::endl;
		TS_ASSERT_DELTA( -22.49034811587382, end_score, 1e-12 );
	}

	void test_etrable_trie_vs_trie() {
		using namespace core::graph;
		using namespace core::pose;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;
		using namespace core::pack;
		using namespace core::pack::rotamer_set;
		using namespace core::pack::task;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn( pose );

		EnergyMethodOptions options; // default is fine
		TableLookupEtableEnergy etab_energy( *( ScoringManager::get_instance()->etable( options.etable_type() ).lock()), options );

		PackerTaskOP task = TaskFactory::create_packer_task( pose );
		for ( Size ii = 1; ii <= 7; ++ii ) task->nonconst_residue_task(ii).prevent_repacking();
		for ( Size ii = 12; ii <= pose.size(); ++ii ) task->nonconst_residue_task(ii).prevent_repacking();
		for ( Size ii = 8; ii <= 11; ++ii ) task->nonconst_residue_task( ii ).or_ex1( true );
		for ( Size ii = 8; ii <= 11; ++ii ) task->nonconst_residue_task( ii ).or_ex2( true );

		GraphOP packer_neighbor_graph = create_packer_graph( pose, sfxn, task );
		RotamerSets rotsets; rotsets.set_task( task );
		rotsets.build_rotamers( pose, sfxn, packer_neighbor_graph );

		// Create tries for each of the three rotamer sets
		for ( Size ii = 8; ii <= 11; ++ii ) etab_energy.prepare_rotamers_for_packing( pose, *rotsets.rotamer_set_for_residue( ii ) );
		for ( Size ii = 8; ii <= 11; ++ii ) TS_ASSERT( rotsets.rotamer_set_for_residue(ii)->get_trie( etable_method ).get() != 0 );

		Size count_comparisons( 0 );
		for ( Size ii = 8; ii <= 11; ++ii ) {
			RotamerSet const & iiset = *rotsets.rotamer_set_for_residue( ii );
			for ( Size jj = ii+1; jj <= 11; ++jj ) {
				RotamerSet const & jjset = *rotsets.rotamer_set_for_residue( jj );
				// compute the rotamer pair energies for ii/jj
				ObjexxFCL::FArray2D< core::PackerEnergy > energy_table( jjset.num_rotamers(), iiset.num_rotamers(), core::PackerEnergy(0.0) );
				etab_energy.evaluate_rotamer_pair_energies( iiset, jjset, pose, sfxn, sfxn.weights(), energy_table );

				/// And now verify that it worked
				ObjexxFCL::FArray2D< core::PackerEnergy > temp_table3( energy_table );
				temp_table3 = 0;
				EnergyMap emap;
				for ( Size kk = 1, kk_end = iiset.num_rotamers(); kk <= kk_end; ++kk ) {
					for ( Size ll = 1, ll_end = jjset.num_rotamers(); ll <= ll_end; ++ll ) {
						++count_comparisons;
						emap.zero();
						etab_energy.residue_pair_energy( *iiset.rotamer( kk ), *jjset.rotamer( ll ), pose, sfxn, emap );
						temp_table3( ll, kk ) += sfxn.weights().dot( emap );
						TS_ASSERT( std::abs( energy_table( ll, kk ) - temp_table3( ll, kk )) <= 1e-3 );
						if ( std::abs( energy_table( ll, kk ) - temp_table3( ll, kk )) > 1e-3 ) {
							std::cout << "Residues " << iiset.resid() << " & " << jjset.resid() << " rotamers: " << kk << " & " << ll;
							std::cout << " tvt/reg discrepancy: tvt= " <<  energy_table( ll, kk ) << " reg= " << temp_table3( ll, kk );
							std::cout << " delta: " << energy_table( ll, kk ) - temp_table3( ll, kk ) << std::endl;
						}
					}
				}
			}
		}
		//std::cout << "Total energy comparisons: " << count_comparisons << std::endl;
	}

	void test_retrieve_etables_from_scoring_manager() {
		methods::EnergyMethodOptions etab_opts1;
		methods::EnergyMethodOptions etab_opts2;
		EtableCAP etable1 = ScoringManager::get_instance()->etable( etab_opts1 );
		EtableCAP etable2 = ScoringManager::get_instance()->etable( etab_opts2 );

		EtableCOP etable1op = etable1.lock();
		EtableCOP etable2op = etable2.lock();

		TS_ASSERT_EQUALS( etable1op, etable2op );
	}

	void test_serialize_etable_trie() {
		TS_ASSERT( true );
#ifdef    SERIALIZATION
		using namespace core::graph;
		using namespace core::pose;
		using namespace core::scoring::etable;
		using namespace core::scoring::etable::etrie;
		using namespace core::scoring::methods;
		using namespace core::scoring::trie;
		using namespace core::pack;
		using namespace core::pack::rotamer_set;
		using namespace core::pack::task;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn( pose );

		PackerTaskOP task = TaskFactory::create_packer_task( pose );
		for ( Size ii = 1; ii <= pose.size(); ++ii ) {
			if ( ii == 8 ) continue;
			task->nonconst_residue_task(ii).prevent_repacking();
		}
		task->nonconst_residue_task( 8 ).or_ex1( true );
		task->nonconst_residue_task( 8 ).or_ex2( true );

		GraphOP packer_neighbor_graph = create_packer_graph( pose, sfxn, task );
		RotamerSets rotsets; rotsets.set_task( task );
		rotsets.build_rotamers( pose, sfxn, packer_neighbor_graph );

		EnergyMethodOptions options; // default is fine
		AnalyticEtableEnergy etab_energy( *( ScoringManager::get_instance()->etable( options.etable_type() ).lock()), options, true );
		RotamerSetOP rotset8 = rotsets.rotamer_set_for_residue( 8 );

		etab_energy.prepare_rotamers_for_packing( pose, *rotset8 );

		conformation::AbstractRotamerTrieCOP trie8 = rotset8->get_trie( methods::etable_method );

		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arc( oss );
			arc( trie8 );
		}

		conformation::AbstractRotamerTrieOP trie8copy;
		std::istringstream iss( oss.str() );
		{
			cereal::BinaryInputArchive arc( iss );
			arc( trie8copy );
		}

		typedef RotamerTrie< EtableAtom, CountPairData_1_2 > etable_rotamer_trie_12;
		typedef std::shared_ptr< etable_rotamer_trie_12 const > etable_rotamer_trie_12_COP;
		etable_rotamer_trie_12_COP orig_trie = utility::pointer::dynamic_pointer_cast< etable_rotamer_trie_12 const > ( trie8 );
		etable_rotamer_trie_12_COP copy_trie = utility::pointer::dynamic_pointer_cast< etable_rotamer_trie_12 const > ( trie8copy );
		TS_ASSERT( trie8copy );

		for ( Size ii = 1; ii <= orig_trie->atoms().size(); ++ii ) {
			TS_ASSERT_EQUALS( orig_trie->atoms()[ ii ].atom(), copy_trie->atoms()[ ii ].atom() );
			TS_ASSERT_EQUALS( orig_trie->atoms()[ ii ].cp_data(), copy_trie->atoms()[ ii ].cp_data() );
			TS_ASSERT_EQUALS( orig_trie->atoms()[ ii ].sibling(), copy_trie->atoms()[ ii ].sibling() );
			TS_ASSERT_EQUALS( orig_trie->atoms()[ ii ].subtree_interaction_sphere_square_radius(), copy_trie->atoms()[ ii ].subtree_interaction_sphere_square_radius() );
			TS_ASSERT_EQUALS( orig_trie->atoms()[ ii ].num_rotamers_in_subtree(), copy_trie->atoms()[ ii ].num_rotamers_in_subtree() );
			TS_ASSERT_EQUALS( orig_trie->atoms()[ ii ].first_atom_in_branch(), copy_trie->atoms()[ ii ].first_atom_in_branch() );
			TS_ASSERT_EQUALS( orig_trie->atoms()[ ii ].has_sibling(), copy_trie->atoms()[ ii ].has_sibling() );
			TS_ASSERT_EQUALS( orig_trie->atoms()[ ii ].is_hydrogen(), copy_trie->atoms()[ ii ].is_hydrogen() );
			TS_ASSERT_EQUALS( orig_trie->atoms()[ ii ].is_rotamer_terminal(), copy_trie->atoms()[ ii ].is_rotamer_terminal() );
		}


#endif // SERIALIZATION

	}

};
