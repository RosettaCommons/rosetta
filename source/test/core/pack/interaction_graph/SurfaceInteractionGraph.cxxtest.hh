// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/SurfaceInteractionGraph.cxxtest.hh
/// @brief  test suite for the surface interaction graph
/// @author Ron Jacak

// Test framework headers
#include <cxxtest/TestSuite.h>

// Core Headers
#include <core/pack/interaction_graph/PDInteractionGraph.hh>
#include <core/pack/interaction_graph/LinearMemoryInteractionGraph.hh>
#include <core/pack/interaction_graph/SurfaceInteractionGraph.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/pack/packer_neighbors.hh>

#include <core/pack/annealer/AnnealerFactory.hh>
#include <core/pack/annealer/SimAnnealerBase.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>

#include <basic/options/keys/packing.OptionKeys.gen.hh>


#include <core/types.hh>


// Auto-header: duplicate removed #include <core/types.hh>

// Utility Headers
#include <basic/Tracer.hh>

// Numeric headers

// Test headers
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

//Auto Headers
#include <core/conformation/AbstractRotamerTrie.fwd.hh>
#include <core/conformation/PointGraph.fwd.hh>
#include <core/conformation/PointGraphData.fwd.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/graph/Graph.hh>
#include <core/graph/UpperEdgeGraph.fwd.hh>
#include <core/graph/unordered_object_pool.fwd.hpp>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/file_data.fwd.hh>
#include <core/pack/annealer/AnnealerFactory.fwd.hh>
#include <core/pack/annealer/SimAnnealerBase.fwd.hh>
#include <core/pack/interaction_graph/AdditionalBackgroundNodesInteractionGraph.hh>
#include <core/pack/interaction_graph/AminoAcidNeighborSparseMatrix.hh>
#include <core/pack/interaction_graph/FixedBBInteractionGraph.fwd.hh>
#include <core/pack/interaction_graph/FixedBBInteractionGraph.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/interaction_graph/LinearMemoryInteractionGraph.fwd.hh>
#include <core/pack/interaction_graph/OnTheFlyInteractionGraph.fwd.hh>
#include <core/pack/interaction_graph/OnTheFlyInteractionGraph.hh>
#include <core/pack/interaction_graph/PDInteractionGraph.fwd.hh>
#include <core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.fwd.hh>
#include <core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.hh>
#include <core/pack/interaction_graph/SparseMatrixIndex.hh>
#include <core/pack/interaction_graph/SurfaceInteractionGraph.fwd.hh>
#include <core/pack/interaction_graph/SurfacePotential.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.fwd.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.hh>
#include <core/pack/rotamer_set/RotamerCouplings.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetsBase.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetsBase.hh>
#include <core/pack/task/IGEdgeReweightContainer.fwd.hh>
#include <core/pack/task/ResfileReader.fwd.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/scoring/ContextGraph.fwd.hh>
#include <core/scoring/ContextGraph.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/NeighborList.fwd.hh>
#include <core/scoring/TenANeighborGraph.fwd.hh>
#include <core/scoring/TwelveANeighborGraph.fwd.hh>
#include <core/scoring/methods/Methods.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1.io.hh>
#include <ObjexxFCL/FArray1A.fwd.hh>
#include <ObjexxFCL/FArray1A.hh>
#include <ObjexxFCL/FArray1P.fwd.hh>
#include <ObjexxFCL/FArray1P.hh>
#include <ObjexxFCL/FArray3.fwd.hh>
#include <ObjexxFCL/FArray3.hh>
#include <ObjexxFCL/FArray3D.fwd.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/StaticIndexRange.fwd.hh>
#include <ObjexxFCL/StaticIndexRange.hh>
#include <assert.h>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/DataCache.fwd.hh>
#include <basic/datacache/DataCache.hh>
#include <basic/options/option.hh>
#include <boost/pool/poolfwd.hpp>


static basic::Tracer TR("test.core.pack.interactiongraph.sig");

using namespace core;
using namespace core::pack;
using namespace core::scoring;


// --------------- Test Class --------------- //

class SurfaceInteractionGraphTests : public CxxTest::TestSuite {

public:

	bool suite_initialized;

	// Shared data elements go here.
	pose::Pose pose;
	rotamer_set::RotamerSetsOP rotsets;
	scoring::ScoreFunctionOP scorefxn;
	graph::GraphOP packer_neighbor_graph;

	// the IG pointers def have to be "global" scope
	interaction_graph::LinearMemorySurfaceInteractionGraphOP lmsolig;
	interaction_graph::PDSurfaceInteractionGraphOP pdsig;
	task::PackerTaskOP designtask;
	annealer::SimAnnealerBaseOP annealer;

	ObjexxFCL::FArray1D_int bestrotamer_at_seqpos;
	PackerEnergy bestenergy;
	PackerEnergy currentenergy, previous_energy_for_node, delta_energy;
	float threshold_for_deltaE_inaccuracy;  // has to be a raw float; otherwise lots of function sigs would have to change


	// --------------- Suite-level Fixture --------------- //

	SurfaceInteractionGraphTests() {
		suite_initialized = false;
	}

	virtual ~SurfaceInteractionGraphTests() {}

	static SurfaceInteractionGraphTests *createSuite() {
		return new SurfaceInteractionGraphTests();
	}

	static void destroySuite( SurfaceInteractionGraphTests *suite ) {
		delete suite;
	}

	void initialize_suite() {

		if ( suite_initialized ) return;
		suite_initialized = true;

		// if the tests are run manually (or one suite at a time), that doesn't mute all of the tracer output by default.  Place
		// a mute here because the interaction graphs generate tons of debugging output (in DEBUG mode anyway).
		core_init_with_additional_options( "-no_optH -mute core.io core.init core.scoring core.mm -restore_pre_talaris_2013_behavior -override_rsd_type_limit" );


		// To create a Surface Interaction Graph object, we need to create a few other objects like a Pose, a ScoreFunction,
		// a PackerTask and a RotamerSets object.  Create all of these objects here in the suite-level fixture since they'll
		// get reused throughout the suite.


		// --- Pose ---
		// since this is a test suite, we don't want to read in PDB files from the command line.  just hardcode the tests to use
		// a predefined test PDB file
		//TR << "Reading in pose..." << std::endl;
		core::import_pose::pose_from_pdb( pose, "core/pack/1l2y_renameH.pdb" );

		// --- PackerTask ---
		// create a custom PackerTask, no extra chi, include current, using the surface score and setting the weight
		designtask = task::TaskFactory::create_packer_task( pose );
		task::parse_resfile(pose, *designtask, "core/pack/interaction_graph/resfile");
		designtask->or_include_current( true );

		// --- ScoreFunction ---
		// create a score function using the standard packer weights
		scorefxn = scoring::get_score_function();
		scorefxn->set_weight( scoring::surface, 0.5 );
		(*scorefxn)( pose );
		pose.update_residue_neighbors();

		// calls setup_for_packing on all of the scoring methods being used. (not sure what that call does though.)
		scorefxn->setup_for_packing( pose, designtask->repacking_residues(), designtask->designing_residues() );

		// --- RotamerSets ---
		rotsets = rotamer_set::RotamerSetsOP( new rotamer_set::RotamerSets() );
		rotsets->set_task( designtask ); // sets the moltenres_2_resid and resid_2_moltenres arrays
		//TR << "Building rotamers..." << std::endl;
		packer_neighbor_graph = create_packer_graph( pose, *scorefxn, designtask );
		rotsets->build_rotamers( pose, *scorefxn, packer_neighbor_graph ); // builds the rotamers
		rotsets->prepare_sets_for_packing( pose, *scorefxn );
		//TR << "\tbuilt " << rotsets->nrotamers() << " rotamers at " << rotsets->nmoltenres() << " positions." << std::endl;


		// Most of the tests in this suite need an interaction graph (more specifically, a Pairwise-decomposable IG).
		// Create the PDIG in the suite-fixture and then for the test fixture, just call blanket_assign_state_0.
		// That will "reset" the interaction graph to the clean state for each test. This saves alot of time because
		// the expensive graph creation (including creating nodes and edges and subsequently dropping edges) only happens
		// once.


		// --- InteractionGraph ---
		//TR << "Instantiating PDSurfaceInteractionGraph..." << std::endl;
		pdsig = interaction_graph::PDSurfaceInteractionGraphOP( new interaction_graph::PDSurfaceInteractionGraph( designtask->num_to_be_packed() ) );
		pdsig->set_pose( pose );
		pdsig->set_packer_task( *designtask );
		pdsig->set_rotamer_sets( *rotsets );

		// compute_energies() does some initialization of the interaction graph and computes the energies
		rotsets->compute_energies( pose, *scorefxn, packer_neighbor_graph, static_cast< interaction_graph::InteractionGraphBaseOP >(pdsig) );


		// Now that we have an interaction graph, a pose, scorefunction, etc, we have everything we need to run the
		// packer except for an annealer. Use just a plain FixbbAnnealer. Go ahead and create a FixbbSA here.  In the
		// test case that uses linmem_ig we'll have to recreate the annealer but this state is common to the rest of
		// the tests since they all use a standard PD IG.


		// initialize some other variables that are used in the SimAnnealers constructor
		bestrotamer_at_seqpos.dimension( pose.total_residue() );
		bool start_with_current = false;
		bool calc_rot_freq = false;
		ObjexxFCL::FArray1D_int current_rot_index; current_rot_index.dimension( pose.total_residue(), 0 );
		ObjexxFCL::FArray1D< PackerEnergy > rot_freq; rot_freq.dimension( pdsig->get_num_total_states(), 0.0 );
		utility::vector0<int> rot_to_pack;

		annealer = annealer::AnnealerFactory::create_annealer(
			designtask, rot_to_pack, bestrotamer_at_seqpos, bestenergy, start_with_current, pdsig, rotsets, current_rot_index, calc_rot_freq, rot_freq );

		// temperature isn't so important, but to make things easy use the SA setup_temp() method
		ObjexxFCL::FArray1D_float loopenergy( 500, 0.0 );  // hardcore the number of loops for this array to the maxnumberofouteriterations
		annealer->setup_temperature( loopenergy, 1 );  // 1 would be the first iteration of outer loop
		threshold_for_deltaE_inaccuracy = std::sqrt( annealer->get_temperature() );

	}


	// --------------- Test Fixture --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case.

	void setUp() {
		initialize_suite();
		pdsig->prepare_for_simulated_annealing();
		pdsig->blanket_assign_state_0();
		pdsig->set_errorfull_deltaE_threshold( threshold_for_deltaE_inaccuracy );
	}

	// Shared finalization goes here.
	// All memory allocated via OPs; objects should destroy themselves so nothing else to do here.
	void tearDown() {}


public:

	// --------------- Test Cases --------------- //


	/// @details
	/// Tests to make sure when doing a design on only some residues that certain residues are indeed being treated and set
	/// as background nodes. If this array returns the wrong indices, background nodes are not being set properly.
	///
	void test_bg_node_2_resid() {

		//TR << "Running test_bg_node_2_resid..." << std::endl;
		TS_ASSERT( pdsig->bg_node_2_resid(2) == 2 );
		TS_ASSERT( pdsig->bg_node_2_resid(4) == 9 );
		TS_ASSERT( pdsig->bg_node_2_resid(7) == 13 );
	}


	/// @details
	/// Tests the function consider_substitution() (hereafter cs()).
	/// cs() takes a position and a new state and tells all nb'ing nodes to update their hASA and return a (possibly
	/// inaccurate) estimate of the change in energy.  The energy is incorrect if the IG decideds to procrastinate the
	/// calculation.  This test needs to verify that nodes are updating their hASAs and that the delta energy is
	/// correct.  Test both a mutation that changes a surface node from HP to P and a mutation which does not lead to
	/// any change in the hASA.
	///
	void test_consider_substitution() {

		//TR << "Running test_consider_substitution..." << std::endl;

		// In this new version of the SIG, every substitution causes a recalculation of the surface score. There are no
		// longer "types" of subs.

		// P->P: commit state 382(ARG) on MR 5 (PDB: 8, LYS)
		pdsig->consider_substitution( 5, 382, delta_energy, previous_energy_for_node );
		TS_ASSERT_DELTA( delta_energy, -101.1885, 0.05 );
		pdsig->commit_considered_substitution();

		// HP->P: commit state 186 (TYR) on MR 2 (GLN)
		pdsig->consider_substitution( 2, 186, delta_energy, previous_energy_for_node );
		TS_ASSERT_DELTA( delta_energy, 103.9336, 0.05 );
		pdsig->commit_considered_substitution();

		int states_a[] = { 0, 186, 0, 0, 382, 0 };
		const size_t size = sizeof(states_a) / sizeof (states_a[0]);
		std::vector<int> correct_state_a(states_a, states_a + size);
		TS_ASSERT_EQUALS( pdsig->get_network_state(), correct_state_a );

		// need to also test here that the nodes have updated their hASAs correctly!
		Real hASA_a[] = { 763.9710, 763.9710, 1151.3855, 772.6270, 627.8220, 300.2651, 383.6072, 527.6359, 475.6991, 336.3497 };
		const size_t size_element_a = sizeof(hASA_a) / sizeof(hASA_a[0]);
		std::vector< Real > correct_hASA_vector_a( hASA_a, hASA_a + size_element_a );
		for ( Size ii=0; ii < correct_hASA_vector_a.size(); ++ii ) {
			TS_ASSERT_DELTA( pdsig->get_hASA_for_node_and_nbs( 2 )[ ii ], correct_hASA_vector_a[ ii ], 0.001 );
		}

		// HP(a)->HP(b): commit state 61 PHE on TRP-6 (MR: 3)
		pdsig->consider_substitution( 3, 61, delta_energy, previous_energy_for_node );
		TS_ASSERT_DELTA( delta_energy, -1.3633, 0.05 );
		currentenergy = pdsig->commit_considered_substitution();

		// HP(a)->HP(b): commit state 25 PHE on LEU-7 (MR: 4) // this sub is favorable, but since one position has more than the
		// max allowed patch area, it will come back as an unfavorable sub
		pdsig->consider_substitution( 4, 25, delta_energy, previous_energy_for_node );
		TS_ASSERT_DELTA( delta_energy, -3.5041, 0.05 );
		currentenergy = pdsig->commit_considered_substitution();

		// HP->P: commit state 435 ARG on ILE-4 (MR: 1)
		pdsig->consider_substitution( 1, 435, delta_energy, previous_energy_for_node );
		TS_ASSERT_DELTA( delta_energy, -105.0106, 0.05 );
		currentenergy = pdsig->commit_considered_substitution();

		// need to also test here that the nodes have updated their hASAs correctly!
		Real hASA_b[] = { 676.0080, 676.0080, 1063.4224, 684.6640, 539.8590, 253.1028, 336.4450, 439.6729, 387.7361, 289.1874 };
		const size_t size_element_b = sizeof(hASA_b) / sizeof(hASA_b[0]);
		std::vector< Real > correct_hASA_vector_b( hASA_b, hASA_b + size_element_b );
		for ( Size ii=0; ii < correct_hASA_vector_b.size(); ++ii ) {
			TS_ASSERT_DELTA( pdsig->get_hASA_for_node_and_nbs( 1 )[ ii ], correct_hASA_vector_b[ ii ], 0.001 );
		}
	}


	/// @details
	/// The main things to test with commit sub() is that the total energy returned is correct and that the node counts (for
	/// the changing node *and* all neighboring nodes) are updated.  It's possible that consider doesn't actually compute
	/// the correct energy because of computation procrastination.  That would also mean the node counts would be inaccurate
	/// until after the commit occurred.
	///
	void test_commit_substitution() {

		//TR << "Running test_commit_substitution..." << std::endl;

		// Tests we'll have here:
		// 1) whether the delta_energy returned by the method is correct for the case where the calculation is not
		// procrastinated, and 2) when it is procrastinated.
		// Another thing to test for is that the hASAs reset if one sub is not commit'd but a following sub is.

		// P->HP: commit state 186 (TYR) on MR 2 (GLN)
		pdsig->consider_substitution( 2, 186, delta_energy, previous_energy_for_node );
		TS_ASSERT_DELTA( delta_energy, 5.5484, 0.05 );
		pdsig->commit_considered_substitution();

		int states_a[] = { 0, 186, 0, 0, 0, 0 };
		const size_t size_a = sizeof(states_a) / sizeof (states_a[0]);
		std::vector<int> correct_state_a(states_a, states_a + size_a);
		TS_ASSERT_EQUALS( pdsig->get_network_state(), correct_state_a );

		// need to also test here that the nodes have updated their hASAs correctly!
		Real hASA_a[] = { 794.5879, 794.5879, 1182.0023, 803.2439, 658.4389, 300.2651, 383.6072, 558.2528, 506.3160, 336.3497 };
		const size_t size_element_a = sizeof(hASA_a) / sizeof(hASA_a[0]);
		std::vector< Real > correct_hASA_vector_a( hASA_a, hASA_a + size_element_a );
		for ( Size ii=0; ii < correct_hASA_vector_a.size(); ++ii ) {
			TS_ASSERT_DELTA( pdsig->get_hASA_for_node_and_nbs( 2 )[ ii ], correct_hASA_vector_a[ ii ], 0.001 );
		}


		// HP->P: commit state 435 ARG on ILE-4 (MR: 1)
		pdsig->consider_substitution( 1, 435, delta_energy, previous_energy_for_node );
		// don't commit

		// HP(a)->HP(b): commit state 25 PHE on LEU-7 (MR: 4) // this residue is on the surface so it should have an effect
		pdsig->consider_substitution( 4, 25, delta_energy, previous_energy_for_node );
		TS_ASSERT_DELTA( delta_energy, -4.9673, 0.05 );
		currentenergy = pdsig->commit_considered_substitution();

		int states_b[] = { 0, 186, 0, 25, 0, 0 };
		const size_t size_b = sizeof(states_b) / sizeof (states_b[0]);
		std::vector<int> correct_state_b(states_b, states_b + size_b);
		TS_ASSERT_EQUALS( pdsig->get_network_state(), correct_state_b );


		// --- make a really bad mutation which should procrastinate the surface calculation
		pdsig->set_observed_sufficient_boolean_true(); // this function call ensures that the surface calculation is procrastinated

		pdsig->consider_substitution( 6, 16, delta_energy, previous_energy_for_node );
		// delta_energy should be an estimated value (equal to the PD deltaE alone) if the graph is procrastinating correctly
		TS_ASSERT_DELTA( delta_energy, 192.4929, 5.0 );

		// to make sure the IG procrastinated, check the node state before and after the commit
		int states_c[] = { 0, 186, 0, 25, 0, 0 };
		const size_t size_c = sizeof(states_c) / sizeof (states_c[0]);
		std::vector<int> correct_state_c(states_c, states_c + size_c);
		TS_ASSERT_EQUALS( pdsig->get_network_state(), correct_state_c );

		currentenergy = pdsig->commit_considered_substitution();
		TS_ASSERT_DELTA( currentenergy, 193.2190, 5.0 );

		int states_d[] = { 0, 186, 0, 25, 0, 16 };
		std::vector<int> correct_state_d(states_d, states_d + size_c);
		TS_ASSERT_EQUALS( pdsig->get_network_state(), correct_state_d );

	}

	/// @details
	/// Make some random commits and make sure the total energy current state assignment method is returning the same
	/// thing that was computed in commit_sub().
	///
	void test_get_energy_current_state_assignment() {

		//TR << "Running test_get_energy_current_state_assignment..." << std::endl;

		pdsig->consider_substitution( 2, 106, delta_energy, previous_energy_for_node );
		currentenergy = pdsig->commit_considered_substitution();
		TS_ASSERT( currentenergy == pdsig->get_energy_current_state_assignment() );

		pdsig->consider_substitution( 5, 10, delta_energy, previous_energy_for_node );
		currentenergy = pdsig->commit_considered_substitution();
		TS_ASSERT( currentenergy == pdsig->get_energy_current_state_assignment() );

		pdsig->consider_substitution( 6, 20, delta_energy, previous_energy_for_node );
		currentenergy = pdsig->commit_considered_substitution();
		TS_ASSERT( currentenergy == pdsig->get_energy_current_state_assignment() );

	}

	/// @details
	/// Near the end of sims, lots of rotamers are tried (which change the alt state counts) but then aren't committed.
	/// This test ensures that the graph is resetting state correctly in those cases.
	///
	void test_blanket_reset_alt_state_counts() {

		//TR << "Running test_blanket_reset_alt_state_counts..." << std::endl;

		// need to also test here that the nodes have updated their hASAs correctly!
		Real hASA_before[] = { 740.2496, 1119.0081, 731.5936, 731.5936, 595.4446, 809.9226, 495.2585, 443.3217, 276.4560, 215.9153 };
		const size_t size_element_before = sizeof(hASA_before) / sizeof(hASA_before[0]);
		std::vector< Real > correct_hASA_vector_before( hASA_before, hASA_before + size_element_before );
		for ( Size ii=0; ii < correct_hASA_vector_before.size(); ++ii ) {
			TS_ASSERT_DELTA( pdsig->get_hASA_for_node_and_nbs( 4 )[ ii ], correct_hASA_vector_before[ ii ], 0.001 );
		}

		pdsig->consider_substitution( 4, 1, delta_energy, previous_energy_for_node );

		Real hASA_after[] = { 705.2015, 1083.9600, 696.5455, 696.5455, 560.3965, 774.8744, 460.2104, 408.2736, 241.4079, 180.8671 };
		const size_t size_element_after = sizeof(hASA_after) / sizeof(hASA_after[0]);
		std::vector< Real > correct_hASA_vector_after( hASA_after, hASA_after + size_element_after );
		for ( Size ii=0; ii < correct_hASA_vector_after.size(); ++ii ) {
			TS_ASSERT_DELTA( pdsig->get_alt_state_hASA_for_node_and_nbs( 4 )[ ii ], correct_hASA_vector_after[ ii ], 0.001 );
		}

		// --- consider another P->HP sub at a nearby residue and check the counts
		// the counts should have been reset to what they were before the consider!
		pdsig->consider_substitution( 5, 550, delta_energy, previous_energy_for_node );

		for ( Size ii=0; ii < correct_hASA_vector_before.size(); ++ii ) {
			TS_ASSERT_DELTA( pdsig->get_hASA_for_node_and_nbs( 4 )[ ii ], correct_hASA_vector_before[ ii ], 0.001 );
		}

	}


	/// @brief
	/// a simple packing run that uses surface scoring and a standard PD interaction graph
	///
	void x_test_partial_redesign_using_pd_ig() {

		//TR << "Running test_partial_redesign_using_pd_ig..." << std::endl;

		// Set these really cool options Mike Tyka added to reduce the number of cycles the annealer runs
		// so that this test doesn't take forever.
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		basic::options::option[ basic::options::OptionKeys::packing::outeriterations_scaling ].value( 0.5 );
		basic::options::option[ basic::options::OptionKeys::packing::inneriterations_scaling ].value( 0.2 );

		// --- InteractionGraph ---
		bool start_with_current = false;
		bool calc_rot_freq = false;
		ObjexxFCL::FArray1D_int current_rot_index; current_rot_index.dimension( pose.total_residue(), 0 );
		ObjexxFCL::FArray1D< PackerEnergy > rot_freq; rot_freq.dimension( pdsig->get_num_total_states(), 0.0 );
		utility::vector0<int> rot_to_pack;

		annealer::SimAnnealerBaseOP redesign_annealer = annealer::AnnealerFactory::create_annealer(
			designtask, rot_to_pack, bestrotamer_at_seqpos, bestenergy, start_with_current, pdsig, rotsets, current_rot_index, calc_rot_freq, rot_freq );

		// temperature isn't so important, but to make things easy use the SA setup_temp() method
		ObjexxFCL::FArray1D_float loopenergy( 500, 0.0 );  // hardcore the number of loops for this array to the maxnumberofouteriterations
		annealer->setup_temperature( loopenergy, 1 );  // 1 would be the first iteration of outer loop

		redesign_annealer->run();

		// Don't worry about checking the final energy. I just want the test to run to completion to make sure there
		// are no bugs that show up only over an annealing run.
		//TS_ASSERT( bestenergy < 30.0 );

		/* this section tests scoring, as opposed to testing the interaction graph; thus, leave it out.
		// after the annealer is done runnning, place new rotamers on the input pose
		for ( core::uint ii = 1; ii <= rotsets->nmoltenres(); ++ii ) {
		core::uint iiresid = rotsets->moltenres_2_resid( ii );
		core::uint iibestrot = rotsets->rotid_on_moltenresidue( bestrotamer_at_seqpos( iiresid ) );
		conformation::ResidueCOP bestrot( rotsets->rotamer_set_for_moltenresidue( ii )->rotamer( iibestrot ) );

		conformation::ResidueOP newresidue( bestrot->create_residue() );
		pose.replace_residue( iiresid, *newresidue, false );
		}

		Energy design_score = (*scorefxn)( pose );
		// TS_ASSERT( design_score < native_score ); // do we want a native score?
		TS_ASSERT_DELTA( design_score, 10.0875, 5.0 );
		*/
	}

	/// @brief
	/// A unit (maybe more of an integration test) that does a short design run using the surface score and a linear memory
	/// interaction graph.
	///
	void x_test_partial_redesign_using_linmem_ig() {

		//TR << "Running test_partial_redesign_using_linmem_ig..." << std::endl;

		// Set these really cool options Mike Tyka added to reduce the number of cycles the annealer runs
		// so that this test doesn't take forever.
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		basic::options::option[ basic::options::OptionKeys::packing::outeriterations_scaling ].value( 0.5 );
		basic::options::option[ basic::options::OptionKeys::packing::inneriterations_scaling ].value( 0.2 );

		// setting an option in the code, very bad coding practice. there's no way to set the value of the recent
		// history size via a PackerTask. only way that exists now is to use the command line which we can't do for
		// just this one test and not the other tests in this Suite. so manually set the option here. it's just a unit test!
		basic::options::option[ basic::options::OptionKeys::packing::linmem_ig ].value( 10 );

		// for this test, we don't want to use the PDIG that's created in the fixture

		// TR << "Instantiating LinearMemorySurfaceInteractionGraph..." << std::endl;
		lmsolig = interaction_graph::LinearMemorySurfaceInteractionGraphOP( new interaction_graph::LinearMemorySurfaceInteractionGraph( designtask->num_to_be_packed() ) );
		lmsolig->set_pose( pose );
		lmsolig->set_packer_task( *designtask );
		lmsolig->set_score_function( *scorefxn );
		lmsolig->set_rotamer_sets( *rotsets );

		// lmsolig is a LinearMem SIG OP; compute_energies() wants an IGBase OP. Need to dereference the
		// OP to get the SIG and then cast to an IGBase.

		// compute_energies() does some initialization of the interaction graph and computes the energies
		rotsets->compute_energies( pose, *scorefxn, packer_neighbor_graph, static_cast< interaction_graph::InteractionGraphBaseOP >(lmsolig) );

		/// Parameters passed by reference in annealers constructor to which it writes at the completion of sim annealing.
		bool start_with_current = false;
		bool calc_rot_freq = false;
		ObjexxFCL::FArray1D_int current_rot_index; current_rot_index.dimension( pose.total_residue(), 0 );
		ObjexxFCL::FArray1D< PackerEnergy > linmem_ig_test_rot_freq( lmsolig->get_num_total_states(), 0.0 );
		utility::vector0<int> rot_to_pack;

		annealer::SimAnnealerBaseOP redesign_annealer = annealer::AnnealerFactory::create_annealer(
			designtask, rot_to_pack, bestrotamer_at_seqpos, bestenergy, start_with_current, lmsolig,
			rotsets, current_rot_index, calc_rot_freq, linmem_ig_test_rot_freq );

		redesign_annealer->run();

		// Don't worry about checking the final energy. I just want the test to run to completion to make sure there
		// are no bugs that show up only over an annealing run.
		//TS_ASSERT( bestenergy < 40.0 );

		/* this section tests scoring, as opposed to testing the interaction graph; thus, leave it out.
		// after the annealer is done runnning, place new rotamers on the input pose
		for ( core::uint ii = 1; ii <= rotsets->nmoltenres(); ++ii ) {
		core::uint iiresid = rotsets->moltenres_2_resid( ii );
		core::uint iibestrot = rotsets->rotid_on_moltenresidue( bestrotamer_at_seqpos( iiresid ) );
		conformation::ResidueCOP bestrot( rotsets->rotamer_set_for_moltenresidue( ii )->rotamer( iibestrot ) );

		conformation::ResidueOP newresidue( bestrot->create_residue() );
		pose.replace_residue( iiresid, *newresidue, false );
		}

		Energy design_score = (*scorefxn)( pose );
		TS_ASSERT( design_score < native_score );  // native score?
		TS_ASSERT_DELTA( design_score, 41.3, 5.0 );
		*/

	}

};
