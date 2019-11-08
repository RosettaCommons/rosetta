// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/FASTERInteractionGraph.cxxtest.hh
/// @brief  test suite for the FASTER interaction graph
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org) -- Added unit test for repeated multi-threaded computation of FASTER interaction graph.

// Test framework headers
#include <cxxtest/TestSuite.h>

// Core Headers
#include <core/pack/interaction_graph/FASTERInteractionGraph.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>

#include <utility/graph/Graph.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>


// Test headers
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

#define TEST_DELTA 1e-5

static basic::Tracer TR("test.core.pack.interaction_graph.FASTERInteractionGraph");

class FASTERInteractionGraphTests : public CxxTest::TestSuite {
public:

	void setUp() {
#ifdef MULTI_THREADED
		core_init_with_additional_options( "-multithreading:total_threads 24 -override_rsd_type_limit" );
#else
		core_init_with_additional_options( "-override_rsd_type_limit" );
#endif
	}

	bool faster_graphs_equal(
		core::pack::interaction_graph::FASTERInteractionGraph const & graph1,
		core::pack::interaction_graph::FASTERInteractionGraph const & graph2,
		core::Size const attempt_number
	) {
		bool returnval(true);

		std::string const errmsg( "\tAttempt " + std::to_string(attempt_number) + ": " );
		//Compare node counts:
		core::Size const nodecount( graph1.get_num_nodes() );
		if ( static_cast< core::Size >( graph2.get_num_nodes() ) != nodecount ) {
			TR << errmsg << "Number of nodes don't match!" << std::endl;
			return false;
		}

		//Compare nodes:
		for ( core::Size i(1); i<=nodecount; ++i ) {
			core::Size const numstates( graph1.get_num_states_for_node(i) );
			if ( numstates != static_cast< core::Size >( graph2.get_num_states_for_node(i) ) ) {
				TR << errmsg << "Number of states for node " << i << " don't match!" << std::endl;
				return false;
			}
			for ( core::Size j(1); j<=numstates; ++j ) {
				//I know this is horrible, but I'm going to cast away the constness JUST for the purpose of this unit test.  NOVICE DEVELOPERS: NEVER, NEVER, NEVER DO THIS!!! --VKM, 5 November 2019.
				if ( std::abs( const_cast< core::pack::interaction_graph::FASTERInteractionGraph & >( graph1 ).get_one_body_energy_for_node_state( i, j ) - const_cast< core::pack::interaction_graph::FASTERInteractionGraph & >( graph2 ).get_one_body_energy_for_node_state( i, j ) ) > TEST_DELTA ) {
					TR << errmsg << "Energy for node " << i << ", state " << j << " doesn't match!" << std::endl;
					returnval = false;
				}
			}
		}

		//Compare edges:
		for ( core::Size ii(1); ii<nodecount; ++ii ) {
			for ( core::Size jj(ii+1); jj<=nodecount; ++jj ) {
				core::pack::interaction_graph::FASTEREdge const * const edge1( graph1.get_faster_edge( ii, jj ) );
				core::pack::interaction_graph::FASTEREdge const * const edge2( graph2.get_faster_edge( ii, jj ) );
				if ( edge1 == nullptr && edge2 == nullptr ) continue; //No edge here.

				if ( ( edge1 == nullptr || edge2 == nullptr  ) && ( edge1 != nullptr || edge2 != nullptr ) ) { //Edge missing in one graph only.
					TR << errmsg << "One graph has an edge between nodes " << ii << " and " << jj << ", while the other does not!" << std::endl;
					return false;
				}

				//If we reach here, both graphs have the edge.

				//Loop through all pairs of states and check all twobody energies:
				core::Size const numstates1( graph1.get_num_states_for_node(ii) );
				core::Size const numstates2( graph2.get_num_states_for_node(jj) );
				for ( core::Size iii(1); iii <= numstates1; ++iii ) {
					for ( core::Size jjj(1); jjj <= numstates2; ++jjj ) {
						if ( std::abs( edge1->get_two_body_energy(iii, jjj) - edge2->get_two_body_energy(iii, jjj) ) > TEST_DELTA ) {
							TR << errmsg + "The twobody energy between node " << ii << ", state " << iii << " and node " << jj << ", state " << jjj << " isn't equal between the two graphs!" << std::endl;
							returnval = false;
						}
					}
				}
			}
		}

		return returnval;
	}

	void test_repeated_multithreaded_FASTER_ig_calculation() {
		using namespace core::chemical;
		using namespace utility::graph;
		using namespace core::pack;
		using namespace core::pack::interaction_graph;
		using namespace core::pack::rotamer_set;
		using namespace core::pack::task;
		using namespace core::pose;
		using namespace core::scoring;

#ifdef MULTI_THREADED
		TR << "Starting test_repeated_multithreaded_FASTER_ig_calculation unit test." << std::endl;

		PoseOP trpcage = create_trpcage_ideal_poseop();
		PackerTaskOP task = TaskFactory::create_packer_task( *trpcage );
		utility::vector1< bool > allowed_aas( num_canonical_aas, false );
		allowed_aas[ aa_ala ] = allowed_aas[ aa_phe ] = allowed_aas[ aa_arg ] = allowed_aas[ aa_lys ] = true;
		for ( core::Size ii = 1; ii <= 20; ++ii ) {
			if ( ii == 11 || ii == 12 || ii == 13 || ii == 7 || ii == 3 || ii == 6 ) {
				task->nonconst_residue_task( ii ).restrict_absent_canonical_aas( allowed_aas );
			} else {
				task->nonconst_residue_task( ii ).prevent_repacking();
			}
		}
		ScoreFunctionOP sfxn = core::scoring::get_score_function();
		(*sfxn)( *trpcage ); // score the pose first;
		sfxn->setup_for_packing( *trpcage, task->repacking_residues(), task->designing_residues() );
		GraphOP packer_neighbor_graph = create_packer_graph( *trpcage, *sfxn, task );
		RotamerSetsOP rotsets( utility::pointer::make_shared< RotamerSets >() );
		rotsets->set_task( task );
		rotsets->build_rotamers( *trpcage, *sfxn, packer_neighbor_graph );
		rotsets->prepare_sets_for_packing( *trpcage, *sfxn );

		//First, compute the graph in a single thread:
		FASTERInteractionGraphOP faster_ig_singlethread( utility::pointer::make_shared< FASTERInteractionGraph >( 6 ) );
		rotsets->compute_energies( *trpcage, *sfxn, packer_neighbor_graph, faster_ig_singlethread, 1 );

		core::Size failure_counts(0);

		for ( core::Size i(1); i<=100; ++i ) {
			FASTERInteractionGraphOP faster_ig_multithread( utility::pointer::make_shared< FASTERInteractionGraph >( 6 ) );
			rotsets->compute_energies( *trpcage, *sfxn, packer_neighbor_graph, faster_ig_multithread, 24 );

			if ( !faster_graphs_equal(*faster_ig_singlethread, *faster_ig_multithread, i) ) {
				TR << "Attempt " << i << ": Multithreaded computation FAILED!" << std::endl;
				++failure_counts;
			} else {
				TR << "Attempt " << i << ": Multithreaded computation passed." << std::endl;
			}
		}

		if ( failure_counts == 0 ) {
			TR << "The test_repeated_multithreaded_FASTER_ig_calculation unit test encountered no failures.  This probably means that the FASTER ig is pretty threadsafe.  Hurrah!" << std::endl;
		}
		TS_ASSERT_EQUALS( failure_counts, 0 );

		TR << "Completed test_repeated_multithreaded_FASTER_ig_calculation unit test." << std::endl;
#else
		TR << "The test_repeated_multithreaded_FASTER_ig_calculation unit test does nothing in the single-threaded Rosetta build.  Returning PASS." << std::endl;
#endif
	}

	void test_instantiate_FASTER_ig() {
		using namespace core::chemical;
		using namespace utility::graph;
		using namespace core::pack;
		using namespace core::pack::interaction_graph;
		using namespace core::pack::rotamer_set;
		using namespace core::pack::task;
		using namespace core::pose;
		using namespace core::scoring;
		using core::Size;

		PoseOP trpcage = create_trpcage_ideal_poseop();
		PackerTaskOP task = TaskFactory::create_packer_task( *trpcage );

		utility::vector1< bool > allowed_aas( num_canonical_aas, false );
		allowed_aas[ aa_ala ] = allowed_aas[ aa_phe ] = allowed_aas[ aa_arg ] = true;

		for ( Size ii = 1; ii <= 20; ++ii ) {
			if ( ii == 11 || ii == 12 || ii == 13 ) {
				task->nonconst_residue_task( ii ).restrict_absent_canonical_aas( allowed_aas );
			} else {
				task->nonconst_residue_task( ii ).prevent_repacking();
			}
		}

		ScoreFunctionOP sfxn = core::scoring::get_score_function();
		(*sfxn)( *trpcage ); // score the pose first;
		sfxn->setup_for_packing( *trpcage, task->repacking_residues(), task->designing_residues() );
		GraphOP packer_neighbor_graph = create_packer_graph( *trpcage, *sfxn, task );

		RotamerSetsOP rotsets( new RotamerSets() );
		rotsets->set_task( task );
		rotsets->build_rotamers( *trpcage, *sfxn, packer_neighbor_graph );
		rotsets->prepare_sets_for_packing( *trpcage, *sfxn );


		FASTERInteractionGraphOP faster_ig( new FASTERInteractionGraph( 3 ) );
		//core::pack::pack_rotamers_setup( *trpcage, *sfxn, task, rot_sets, ig );

#ifdef MULTI_THREADED
		rotsets->compute_energies( *trpcage, *sfxn, packer_neighbor_graph, faster_ig, 24 );
#else
		rotsets->compute_energies( *trpcage, *sfxn, packer_neighbor_graph, faster_ig, 1 );
#endif

		/*for ( Size ii = 1; ii <= rotsets->nmoltenres(); ++ii ) {
		std::cout << "Rotset " << ii << " with " << rotsets->rotamer_set_for_moltenresidue(ii)->num_rotamers() << " rotamers" << std::endl;
		}*/

		//std::cout.precision( 8 );

		faster_ig->prepare_graph_for_simulated_annealing();
		faster_ig->prepare_for_FASTER();

		faster_ig->assign_BMEC();
		//faster_ig->print_vertices(); // DETERMINE the BMEC by printing out the one body energies with this function.
		ObjexxFCL::FArray1D_int netstate( 3 );
		faster_ig->get_current_network_state( netstate );
		/*std::cout << "BMEC state:";
		for ( Size ii = 1; ii <= 3; ++ii ) {
		std::cout << " " << netstate( ii );
		}
		std::cout << std::endl;
		std::cout << faster_ig->get_energy_current_state_assignment() << std::endl;*/

		TR << "netstates 1, 2, and 3: " << netstate(1) << ", " << netstate(2) << ", and " <<netstate(3) << std::endl;

		TS_ASSERT( netstate( 1 ) ==  1 );
		TS_ASSERT( netstate( 2 ) ==  1 );
		TS_ASSERT( netstate( 3 ) ==  1 );

		//for ( core::Size ii = 1; ii <= 3; ++ii ) {
		// trpcage->replace_residue( ii + 10, *rotsets->rotamer_set_for_moltenresidue( ii )->rotamer( netstate(ii) ), false );
		//}
		//std::cout << "score for bmec" << (*sfxn)( *trpcage );

		TR << "faster_ig->get_energy_current_state_assignment() = " << TR.precision(10) << faster_ig->get_energy_current_state_assignment() << std::endl;
		TS_ASSERT_DELTA( faster_ig->get_energy_current_state_assignment(), 2.667302132, 1e-5 );

		faster_ig->relax_in_current_context();
		faster_ig->commit_relaxation();
		faster_ig->get_current_network_state( netstate );

		/*std::cout << "BMEC relaxed relaxed state:";
		for ( Size ii = 1; ii <= 3; ++ii ) {
		std::cout << " " << netstate( ii );
		}
		std::cout << std::endl;
		std::cout << faster_ig->get_energy_current_state_assignment() << std::endl;*/

		TR << "netstates 1, 2, and 3: " << netstate(1) << ", " << netstate(2) << ", and " <<netstate(3) << std::endl;
		TS_ASSERT( netstate( 1 ) ==  1 );
		TS_ASSERT( netstate( 2 ) ==  1 );
		TS_ASSERT( netstate( 3 ) ==  1 );

		//std::cout << faster_ig->get_energy_current_state_assignment() << std::endl;

		TR << "faster_ig->get_energy_current_state_assignment() = " << faster_ig->get_energy_current_state_assignment() << std::endl;
		TS_ASSERT_DELTA( faster_ig->get_energy_current_state_assignment(), 2.667302132, 1e-5 );

		faster_ig->perturb_sBR_and_relax( 2, 6 );
		faster_ig->commit_relaxation();
		faster_ig->get_current_network_state( netstate );

		/*std::cout << "sPBR (2,6) relaxed state:";
		for ( Size ii = 1; ii <= 3; ++ii ) {
		std::cout << " " << netstate( ii );
		}
		std::cout << std::endl;
		std::cout << faster_ig->get_energy_current_state_assignment() << std::endl;*/


		TR << "netstates 1, 2, and 3: " << netstate(1) << ", " << netstate(2) << ", and " <<netstate(3) << std::endl;
		TS_ASSERT( netstate( 1 ) ==  1 );
		TS_ASSERT( netstate( 2 ) ==  6 );
		TS_ASSERT( netstate( 3 ) ==  1 );

		TR << "faster_ig->get_energy_current_state_assignment() = " << faster_ig->get_energy_current_state_assignment() << std::endl;
		TS_ASSERT_DELTA( faster_ig->get_energy_current_state_assignment(), 7.454839706, 1e-5 );

		core::PackerEnergy delta1 = faster_ig->perturb_sBR_and_relax( 3, 5 );
		faster_ig->commit_relaxation();
		faster_ig->get_current_network_state( netstate );

		//std::cout << faster_ig->get_energy_current_state_assignment() << std::endl;

		TR << "faster_ig->get_energy_current_state_assignment() = " << faster_ig->get_energy_current_state_assignment() << std::endl;
		TR << "delta1 = " << delta1 << std::endl;
		TS_ASSERT_DELTA( faster_ig->get_energy_current_state_assignment(), 4.343521595, 1e-5 );
		TS_ASSERT_DELTA( delta1, -3.111317396, 1e-5 );

		/*std::cout << "sPBR (3,5) relaxed state:";
		for ( Size ii = 1; ii <= 3; ++ii ) {
		std::cout << " " << netstate( ii );
		}
		std::cout << std::endl;
		std::cout << faster_ig->get_energy_current_state_assignment() << std::endl;*/

		/*core::PackerEnergy delta2 = */faster_ig->perturb_sBR_and_relax( 2, 35 );
		faster_ig->commit_relaxation();
		faster_ig->get_current_network_state( netstate );

		/*std::cout << "sPBR (2,35) relaxed state:";
		for ( Size ii = 1; ii <= 3; ++ii ) {
		std::cout << " " << netstate( ii );
		}
		std::cout << std::endl;
		std::cout << faster_ig->get_energy_current_state_assignment() << std::endl;*/

		TR << "netstates 1, 2, and 3: " << netstate(1) << ", " << netstate(2) << ", and " <<netstate(3) << std::endl;
		TS_ASSERT( netstate( 1 ) ==  1 );
		TS_ASSERT( netstate( 2 ) ==  35 );
		TS_ASSERT( netstate( 3 ) ==  1 );

		TR << "faster_ig->get_energy_current_state_assignment() = " << faster_ig->get_energy_current_state_assignment() << std::endl;
		TS_ASSERT_DELTA( faster_ig->get_energy_current_state_assignment(), 5.77839756, 1e-5 );

		core::PackerEnergy delta3, dummy;
		faster_ig->consider_substitution( 2, 6, delta3, dummy );
		TR << "delta3 = " << delta3 << std::endl;
		TS_ASSERT_DELTA( delta3, 1.676442146, 1e-5 );

		core::PackerEnergy delta4 = faster_ig->perturb_sBR_and_relax( 3, 20 );
		faster_ig->commit_relaxation();
		faster_ig->get_current_network_state( netstate );

		/*std::cout << "sPBR (3,20) relaxed state:";
		for ( Size ii = 1; ii <= 3; ++ii ) {
		std::cout << " " << netstate( ii );
		}
		std::cout << std::endl;
		std::cout << faster_ig->get_energy_current_state_assignment() << std::endl;*/

		TR << "netstates 1, 2, and 3: " << netstate(1) << ", " << netstate(2) << ", and " <<netstate(3) << std::endl;
		TS_ASSERT( netstate( 1 ) ==  1 );
		TS_ASSERT( netstate( 2 ) ==  1 );
		TS_ASSERT( netstate( 3 ) ==  20 );

		TR << "delta4 = " << delta4 << std::endl;
		TS_ASSERT_DELTA( delta4, -0.8559448719, 1e-5 );
		TR << "faster_ig->get_energy_current_state_assignment() = " << faster_ig->get_energy_current_state_assignment() << std::endl;
		TS_ASSERT_DELTA( faster_ig->get_energy_current_state_assignment(), 4.922451973, 1e-5 );


	}


};
