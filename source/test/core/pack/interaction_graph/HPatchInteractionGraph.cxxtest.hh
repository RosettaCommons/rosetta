// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/HPatchInteractionGraph.cxxtest.hh
/// @brief  test suite for the hpatch score based interaction graph
/// @author Ron Jacak

// Test framework headers
#include <cxxtest/TestSuite.h>

// Core Headers
#include <core/conformation/Residue.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/Pose.hh>

#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <core/pack/annealer/AnnealerFactory.hh>
#include <core/pack/annealer/SimAnnealerBase.hh>
#include <core/pack/interaction_graph/PDInteractionGraph.hh>
#include <core/pack/interaction_graph/LinearMemoryInteractionGraph.hh>
#include <core/pack/interaction_graph/HPatchInteractionGraph.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>


#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/sasa.hh>

#include <core/types.hh>

// ObjexxFCL Header

// Utility Headers
#include <basic/Tracer.hh>

// Numeric headers

// Test headers
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

//Auto Headers
#include <core/graph/Graph.hh>
#include <core/import_pose/import_pose.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1.io.hh>
#include <assert.h>
#include <basic/options/option.hh>


static basic::Tracer TR("test.core.pack.interactiongraph.hig");

using namespace core;
using namespace core::pack;
using namespace core::scoring;
using namespace ObjexxFCL::format;


// --------------- Test Class --------------- //

class HPatchInteractionGraphTests : public CxxTest::TestSuite {

public:

	bool suite_initialized;

	// Shared data elements go here.
	pose::Pose pose;
	rotamer_set::RotamerSetsOP rotsets;
	scoring::ScoreFunctionOP scorefxn;
	graph::GraphOP packer_neighbor_graph;

	// the IG pointers def have to be "global" scope
	interaction_graph::LinearMemoryHPatchInteractionGraphOP lmhig;
	interaction_graph::PDHPatchInteractionGraphOP pdhig;
	task::PackerTaskOP designtask;
	annealer::SimAnnealerBaseOP annealer;

	ObjexxFCL::FArray1D_int bestrotamer_at_seqpos;
	PackerEnergy bestenergy;
	PackerEnergy currentenergy, previous_energy_for_node, delta_energy;
	float threshold_for_deltaE_inaccuracy;  // has to be a raw float; otherwise lots of function sigs would have to change


	// --------------- Suite-level Fixture --------------- //

	HPatchInteractionGraphTests() {
		suite_initialized = false;
	}

	virtual ~HPatchInteractionGraphTests() {}

	static HPatchInteractionGraphTests *createSuite() {
		return new HPatchInteractionGraphTests();
	}

	static void destroySuite( HPatchInteractionGraphTests *suite ) {
		delete suite;
	}

	void initialize_suite() {

		if ( suite_initialized ) return;
		suite_initialized = true;

		// if the tests are run manually (or one suite at a time), that doesn't mute all of the tracer output by default.  Place
		// a mute here because the interaction graphs generate tons of debugging output (in DEBUG mode anyway).
		core_init_with_additional_options( "-no_optH -mute core.io core.init core.mm -restore_pre_talaris_2013_behavior -override_rsd_type_limit" );


		// To create a HPatch Interaction Graph object, we need to create a few other objects like a Pose, a ScoreFunction,
		// a PackerTask and a RotamerSets object.  Create all of these objects here in the suite-level fixture since they'll
		// get reused throughout the suite.


		// --- Pose ---
		// since this is a test suite, we don't want to read in PDB files from the command line.  just hardcode the tests to use
		// a predefined test PDB file
		TR << "Reading in pose..." << std::endl;
		core::import_pose::pose_from_pdb( pose, "core/pack/1l2y_renameH.pdb" );

		// --- PackerTask ---
		// create a custom PackerTask, no extra chi, include current, using the surface score and setting the weight
		designtask = task::TaskFactory::create_packer_task( pose );
		task::parse_resfile(pose, *designtask, "core/pack/interaction_graph/resfile");
		designtask->or_include_current( true );

		// --- ScoreFunction ---
		// create a score function using the standard packer weights
		scorefxn = scoring::get_score_function();
		(*scorefxn)( pose );
		pose.update_residue_neighbors();

		// calls setup_for_packing on all of the scoring methods being used. (not sure what that call does though.)
		scorefxn->setup_for_packing( pose, designtask->repacking_residues(), designtask->designing_residues() );

		// --- RotamerSets ---
		rotsets = rotamer_set::RotamerSetsOP( new rotamer_set::RotamerSets() );
		rotsets->set_task( designtask ); // sets the moltenres_2_resid and resid_2_moltenres arrays
		TR << "Building rotamers..." << std::endl;
		packer_neighbor_graph = create_packer_graph( pose, *scorefxn, designtask );
		rotsets->build_rotamers( pose, *scorefxn, packer_neighbor_graph ); // builds the rotamers
		rotsets->prepare_sets_for_packing( pose, *scorefxn );
		TR << "\tbuilt " << rotsets->nrotamers() << " rotamers at " << rotsets->nmoltenres() << " positions." << std::endl;


		// Most of the tests in this suite need an interaction graph (more specifically, a Pairwise-decomposable IG).
		// Create the PDIG in the suite-fixture and then for the test fixture, just call blanket_assign_state_0.
		// That will "reset" the interaction graph to the clean state for each test. This saves alot of time because
		// the expensive graph creation (including creating nodes and edges and subsequently dropping edges) only happens
		// once.


		// --- InteractionGraph ---
		TR << "Instantiating PDHPatchInteractionGraph..." << std::endl;
		pdhig = interaction_graph::PDHPatchInteractionGraphOP( new interaction_graph::PDHPatchInteractionGraph( designtask->num_to_be_packed() ) );
		pdhig->set_pose( pose );
		pdhig->set_packer_task( *designtask );
		pdhig->set_rotamer_sets( *rotsets );

		// compute_energies() does some initialization of the interaction graph and computes the energies
		rotsets->compute_energies( pose, *scorefxn, packer_neighbor_graph, static_cast< interaction_graph::InteractionGraphBaseOP >(pdhig) );


		// Now that we have an interaction graph, a pose, scorefunction, etc, we have everything we need to run the
		// packer except for an annealer. Use just a plain FixbbAnnealer. Go ahead and create a FixbbSA here.  In the
		// test case that uses linmem_ig we'll have to recreate the annealer but this state is common to the rest of
		// the tests since they all use a standard PD IG.


		// initialize some other variables that are used in the SimAnnealers constructor
		bestrotamer_at_seqpos.dimension( pose.total_residue() );
		bool start_with_current = false;
		bool calc_rot_freq = false;
		ObjexxFCL::FArray1D_int current_rot_index; current_rot_index.dimension( pose.total_residue(), 0 );
		ObjexxFCL::FArray1D< PackerEnergy > rot_freq; rot_freq.dimension( pdhig->get_num_total_states(), 0.0 );
		utility::vector0<int> rot_to_pack;

		annealer = annealer::AnnealerFactory::create_annealer(
			designtask, rot_to_pack, bestrotamer_at_seqpos, bestenergy, start_with_current, utility::pointer::dynamic_pointer_cast<interaction_graph::AnnealableGraphBase>(pdhig), rotsets, current_rot_index, calc_rot_freq, rot_freq );

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
		TR << "Initializing IG..." << std::endl;
		pdhig->prepare_for_simulated_annealing();
		pdhig->blanket_assign_state_0();
		pdhig->set_errorfull_deltaE_threshold( threshold_for_deltaE_inaccuracy );
		//std::cout << std::endl;
	}

	// Shared finalization goes here.
	// All memory allocated via OPs; objects should destroy themselves so nothing else to do here.
	void tearDown() {
		//TR << "called destructor" << std::endl;
	}


public:

	// --------------- Test Cases --------------- //


	/// @details
	/// Tests to make sure when doing a design on only some residues that certain residues are indeed being treated and set
	/// as background nodes. If this array returns the wrong indices, background nodes are not being set properly.
	///
	void test_bg_node_2_resid() {
		TR << "Running test_bg_node_2_resid..." << std::endl;
		TS_ASSERT( pdhig->bg_node_2_resid(2) == 2 );
		TS_ASSERT( pdhig->bg_node_2_resid(4) == 9 );
		TS_ASSERT( pdhig->bg_node_2_resid(7) == 13 );
	}

	/// @details
	/// Tests to make sure the graph is properly calculating SASA for all the nodes and bgnodes.
	/// Obtains the sasa for all nodes and bgnodes when calculated with the RotamerDots. Then we compare those values
	/// to what gets computed for the residues by calc_per_atom_sasa in sasa.cc.
	///
	void test_graph_initialization() {
		TR << "Running test_graph_initialization..." << std::endl;
		utility::vector1< Real > node_sasas( (Size)pdhig->get_num_nodes(), 0.0 );
		utility::vector1< Real > bgnode_sasas ( (Size)pdhig->get_num_background_nodes(), 0.0 );
		pdhig->get_all_sasas( node_sasas, bgnode_sasas );

		id::AtomID_Map< Real > atom_sasa;
		utility::vector1< Real > rsd_sasa;
		Real probe_radius = 1.4;

		// create a very custom atom_subset mask
		//
		// manually disable inclusion of residues 4-8 and 11 (which are the residues which *are* designable in this test case).
		// they will be FC Nodes that are in the unassigned state so they produce no sasa overlap.
		// to get the right sasas on the BG Nodes we have to exclude them from the sasa calculations
		//
		id::AtomID_Map< bool > atom_subset;
		atom_subset.resize( pose.n_residue() );
		// init all heavy atoms to true and all H's to false
		for ( Size ii=1; ii <= pose.n_residue(); ++ii ) {
			atom_subset.resize( ii, pose.residue_type(ii).natoms(), false );
			for ( Size jj = 1; jj <= pose.residue_type(ii).nheavyatoms(); ++jj ) {
				atom_subset[ ii ][ jj ] = true;
			}
		}
		for ( Size ii=4; ii <= 8; ++ii ) {
			for ( Size jj=1; jj <= pose.residue_type(ii).nheavyatoms(); ++jj ) {
				atom_subset[ ii ][ jj ] = false;
			}
		}
		for ( Size jj=1; jj <= pose.residue_type(11).nheavyatoms(); ++jj ) {
			atom_subset[ 11 ][ jj ] = false;
		}

		core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius, false /* no big polar H */,
			atom_subset, true /* use_naccess_sasa_radii */, true /* expand polar radii */ );

		// all the Nodes should have 0.0 for their SASA because they should all be in the unassigned state after initialization
		for ( Size ii=1; ii <= node_sasas.size(); ++ii ) {
			TS_ASSERT_DELTA( node_sasas[ ii ], 0.0, 0.01 );
		}
		// all the background nodes should have some value for SASA
		for ( Size ii=1; ii <= bgnode_sasas.size(); ++ii ) {
			TS_ASSERT_DELTA( bgnode_sasas[ ii ], rsd_sasa[ pdhig->bg_node_2_resid( ii ) ], 0.01 );
		}

	}

	/// @details
	/// Tests the function get_curr_state_hphobes() as well as the alternate state one.
	void test_curr_state_hphobes() {
		TR << "Running test_curr_state_hphobes..." << std::endl;

		// P->P: commit state 382(ARG) on MR 5 (PDB: 8, LYS)
		// this sub only produces SASA changes on bg node indexes 4 and 5 (1-3 don't change)
		pdhig->consider_substitution( 5, 382, delta_energy, previous_energy_for_node );
		pdhig->commit_considered_substitution();

		// HP->P: commit state 186 (TYR) on MR 2 (GLN)
		pdhig->consider_substitution( 2, 186, delta_energy, previous_energy_for_node );
		pdhig->commit_considered_substitution();

		// HP(a)->HP(b): commit state 61 PHE on TRP-6 (MR: 3) #ha: 160
		pdhig->consider_substitution( 3, 61, delta_energy, previous_energy_for_node );
		pdhig->commit_considered_substitution();

		// HP(a)->P(b): commit state 25 GLU on LEU-7 (MR: 4) #ha: 161
		pdhig->consider_substitution( 4, 25, delta_energy, previous_energy_for_node );
		pdhig->commit_considered_substitution();

		// HP->P: commit state 435 ARG on ILE-4 (MR: 1) #ha: 164
		pdhig->consider_substitution( 1, 435, delta_energy, previous_energy_for_node );
		pdhig->commit_considered_substitution();

		// commit state 17 SER on GLY-11 (MR: 6)
		pdhig->consider_substitution( 6, 17, delta_energy, previous_energy_for_node );
		pdhig->commit_considered_substitution();

		// all the Nodes are in an assigned state now.

		// the function curr_state_hphobes() returns a vector1 of atom indexes of the hydrophobic atoms in the residue.
		// we can verify the correctness by asking the pdhig to give us the rotamer and checking the C or S status manually.

		for ( Size ii = 1; ii <= (Size)pdhig->get_num_nodes(); ++ii ) {
			utility::vector1< Size > const & ii_node_curr_state_hphobes = pdhig->get_hpatch_node( ii )->curr_state_hphobes();
			conformation::ResidueCOP rotamer = pdhig->get_hpatch_node( ii )->curr_state_rotamer();

			//TR << "test_curr_state_hphobes(): node " << ii << ": ii_node_curr_state_hphobes: [ ";
			for ( Size jj=1; jj <= ii_node_curr_state_hphobes.size(); ++jj ) {
				//TR << ii_node_curr_state_hphobes[ jj ] << ", ";
				TS_ASSERT( ( rotamer->atom_type( ii_node_curr_state_hphobes[ jj ] ).element() == "C" ) ||
					( rotamer->atom_type( ii_node_curr_state_hphobes[ jj ] ).element() == "S" ) );
			}
			//TR << "], ";

			utility::vector1< Size > const & ii_node_alt_state_hphobes = pdhig->get_hpatch_node( ii )->alt_state_hphobes();
			rotamer = pdhig->get_hpatch_node( ii )->alt_state_rotamer();
			//TR << "ii_node_alt_state_hphobes: [ ";
			for ( Size jj=1; jj <= ii_node_alt_state_hphobes.size(); ++jj ) {
				//TR << ii_node_alt_state_hphobes[ jj ] << ", ";
				TS_ASSERT( ( rotamer->atom_type( ii_node_alt_state_hphobes[ jj ] ).element() == "C" ) ||
					( rotamer->atom_type( ii_node_alt_state_hphobes[ jj ] ).element() == "S" ) );
			}
			//TR << "]" << std::endl;
		}


		// commit state 551(TRP) on MR 5 (PDB: 8, LYS)
		pdhig->consider_substitution( 5, 551, delta_energy, previous_energy_for_node );
		//currentenergy = pdhig->commit_considered_substitution();

		for ( Size ii = 1; ii <= (Size)pdhig->get_num_nodes(); ++ii ) {
			utility::vector1< Size > const & ii_node_curr_state_hphobes = pdhig->get_hpatch_node( ii )->curr_state_hphobes();
			conformation::ResidueCOP rotamer = pdhig->get_hpatch_node( ii )->curr_state_rotamer();
			//TR << "test_curr_state_hphobes(): node " << ii << ": ii_node_curr_state_hphobes: [ ";
			for ( Size jj=1; jj <= ii_node_curr_state_hphobes.size(); ++jj ) {
				//TR << ii_node_curr_state_hphobes[ jj ] << ", ";
				TS_ASSERT( ( rotamer->atom_type( ii_node_curr_state_hphobes[ jj ] ).element() == "C" ) ||
					( rotamer->atom_type( ii_node_curr_state_hphobes[ jj ] ).element() == "S" ) );
			}
			//TR << "], ";

			utility::vector1< Size > const & ii_node_alt_state_hphobes = pdhig->get_hpatch_node( ii )->alt_state_hphobes();
			rotamer = pdhig->get_hpatch_node( ii )->alt_state_rotamer();
			//TR << "ii_node_alt_state_hphobes: [ ";
			for ( Size jj=1; jj <= ii_node_alt_state_hphobes.size(); ++jj ) {
				//TR << ii_node_alt_state_hphobes[ jj ] << ", ";
				TS_ASSERT( ( rotamer->atom_type( ii_node_alt_state_hphobes[ jj ] ).element() == "C" ) ||
					( rotamer->atom_type( ii_node_alt_state_hphobes[ jj ] ).element() == "S" ) );
			}
			//TR << "]" << std::endl;
		}

	}


	/// @details
	/// Tests the function get_curr_state_exp_hphobes() as well as the alternate state one. This test is similar to the
	/// one above but this one tests whether the determination of exposed or not is happening correctly.
	///
	void test_curr_state_exp_hphobes() {
		TR << "Running test_curr_state_exp_hphobes..." << std::endl;

		// P->P: commit state 382(ARG) on MR 5 (PDB: 8, LYS)
		// this sub only produces SASA changes on bg node indexes 4 and 5 (1-3 don't change)
		pdhig->consider_substitution( 5, 382, delta_energy, previous_energy_for_node ); pdhig->commit_considered_substitution();

		// HP->P: commit state 186 (TYR) on MR 2 (GLN)
		pdhig->consider_substitution( 2, 186, delta_energy, previous_energy_for_node ); pdhig->commit_considered_substitution();

		// HP(a)->HP(b): commit state 61 PHE on TRP-6 (MR: 3) #ha: 160
		pdhig->consider_substitution( 3, 61, delta_energy, previous_energy_for_node ); pdhig->commit_considered_substitution();

		// HP(a)->P(b): commit state 25 GLU on LEU-7 (MR: 4) #ha: 161
		pdhig->consider_substitution( 4, 25, delta_energy, previous_energy_for_node ); pdhig->commit_considered_substitution();

		// HP->P: commit state 435 ARG on ILE-4 (MR: 1) #ha: 164
		pdhig->consider_substitution( 1, 435, delta_energy, previous_energy_for_node ); pdhig->commit_considered_substitution();

		// commit state 17 SER on GLY-11 (MR: 6)
		pdhig->consider_substitution( 6, 17, delta_energy, previous_energy_for_node ); pdhig->commit_considered_substitution();

		// all the Nodes are in an assigned state now.

		// the function curr_state_exp_hphobes() returns a vector1 of atom indexes of the exposed hydrophobic atoms in the residue.
		// we can verify the correctness by asking the pdhig to give us the rotamer and checking the C or S status manually,
		// and checking that the SASA is > 0.

		// FC nodes, current state and alt state
		for ( Size ii = 1; ii <= (Size)pdhig->get_num_nodes(); ++ii ) {
			utility::vector1< Size > const & ii_node_curr_state_exp_hphobes = pdhig->get_hpatch_node( ii )->curr_state_exp_hphobes();
			conformation::ResidueCOP rotamer = pdhig->get_hpatch_node( ii )->curr_state_rotamer();
			//TR << "test_curr_state_exp_hphobes(): node " << ii << ": ii_node_curr_state_exp_hphobes: [ ";
			for ( Size jj=1; jj <= ii_node_curr_state_exp_hphobes.size(); ++jj ) {
				//TR << ii_node_curr_state_exp_hphobes[ jj ] << ", ";
				TS_ASSERT( ( rotamer->atom_type( ii_node_curr_state_exp_hphobes[ jj ] ).element() == "C" ) ||
					( rotamer->atom_type( ii_node_curr_state_exp_hphobes[ jj ] ).element() == "S" ) );
				TS_ASSERT( pdhig->get_hpatch_node( ii )->get_current_state_sasa( ii_node_curr_state_exp_hphobes[ jj ] ) > 0.0 );
			}
			//TR << "], ";

			utility::vector1< Size > const & ii_node_alt_state_exp_hphobes = pdhig->get_hpatch_node( ii )->alt_state_exp_hphobes();
			rotamer = pdhig->get_hpatch_node( ii )->alt_state_rotamer();
			//TR << "test_curr_state_exp_hphobes(): ii_node_alt_state_exp_hphobes: [ ";
			for ( Size jj=1; jj <= ii_node_alt_state_exp_hphobes.size(); ++jj ) {
				//TR << ii_node_alt_state_exp_hphobes[ jj ] << ", ";
				TS_ASSERT( ( rotamer->atom_type( ii_node_alt_state_exp_hphobes[ jj ] ).element() == "C" ) ||
					( rotamer->atom_type( ii_node_alt_state_exp_hphobes[ jj ] ).element() == "S" ) );
				TS_ASSERT( pdhig->get_hpatch_node( ii )->get_alternate_state_sasa( ii_node_alt_state_exp_hphobes[ jj ] ) > 0.0 );
			}
			//TR << "]" << std::endl;
		}

		// BG nodes, current state and alt state
		for ( Size ii = 1; ii <= (Size)pdhig->get_num_background_nodes(); ++ii ) {
			utility::vector1< Size > const & ii_bg_node_curr_state_exp_hphobes = pdhig->get_hpatch_bg_node( ii )->curr_state_exp_hphobes();
			conformation::ResidueCOP rotamer = pdhig->get_hpatch_bg_node( ii )->get_rotamer();
			//TR << "test_curr_state_exp_hphobes(): bg node " << ii << ": ii_bg_node_curr_state_exp_hphobes: [ ";
			for ( Size jj=1; jj <= ii_bg_node_curr_state_exp_hphobes.size(); ++jj ) {
				//TR << ii_bg_node_curr_state_exp_hphobes[ jj ] << ", ";
				TS_ASSERT( ( rotamer->atom_type( ii_bg_node_curr_state_exp_hphobes[ jj ] ).element() == "C" ) ||
					( rotamer->atom_type( ii_bg_node_curr_state_exp_hphobes[ jj ] ).element() == "S" ) );
				TS_ASSERT( pdhig->get_hpatch_bg_node( ii )->get_current_sasa( ii_bg_node_curr_state_exp_hphobes[ jj ] ) > 0.0 );
			}
			//TR << "], ";

			utility::vector1< Size > const & ii_bg_node_alt_state_exp_hphobes = pdhig->get_hpatch_bg_node( ii )->alt_state_exp_hphobes();
			//TR << "test_curr_state_exp_hphobes(): ii_bg_node_alt_state_exp_hphobes: [ ";
			for ( Size jj=1; jj <= ii_bg_node_alt_state_exp_hphobes.size(); ++jj ) {
				//TR << ii_bg_node_alt_state_exp_hphobes[ jj ] << ", ";
				TS_ASSERT( ( rotamer->atom_type( ii_bg_node_alt_state_exp_hphobes[ jj ] ).element() == "C" ) ||
					( rotamer->atom_type( ii_bg_node_alt_state_exp_hphobes[ jj ] ).element() == "S" ) );
				TS_ASSERT( pdhig->get_hpatch_bg_node( ii )->get_alternate_sasa( ii_bg_node_alt_state_exp_hphobes[ jj ] ) > 0.0 );
			}
			//TR << "]" << std::endl;
		}

		// commit state 551(TRP) on MR 5 (PDB: 8, LYS)
		pdhig->consider_substitution( 5, 551, delta_energy, previous_energy_for_node );

		// FC nodes, current state and alt state
		for ( Size ii = 1; ii <= (Size)pdhig->get_num_nodes(); ++ii ) {
			utility::vector1< Size > const & ii_node_curr_state_exp_hphobes = pdhig->get_hpatch_node( ii )->curr_state_exp_hphobes();
			conformation::ResidueCOP rotamer = pdhig->get_hpatch_node( ii )->curr_state_rotamer();
			//TR << "test_curr_state_exp_hphobes(): node " << ii << ": ii_node_curr_state_exp_hphobes: [ ";
			for ( Size jj=1; jj <= ii_node_curr_state_exp_hphobes.size(); ++jj ) {
				//TR << ii_node_curr_state_exp_hphobes[ jj ] << ", ";
				TS_ASSERT( ( rotamer->atom_type( ii_node_curr_state_exp_hphobes[ jj ] ).element() == "C" ) ||
					( rotamer->atom_type( ii_node_curr_state_exp_hphobes[ jj ] ).element() == "S" ) );
				TS_ASSERT( pdhig->get_hpatch_node( ii )->get_current_state_sasa( ii_node_curr_state_exp_hphobes[ jj ] ) > 0.0 );
			}
			//TR << "], ";

			utility::vector1< Size > const & ii_node_alt_state_exp_hphobes = pdhig->get_hpatch_node( ii )->alt_state_exp_hphobes();
			rotamer = pdhig->get_hpatch_node( ii )->alt_state_rotamer();
			//TR << "test_curr_state_exp_hphobes(): ii_node_alt_state_exp_hphobes: [ ";
			for ( Size jj=1; jj <= ii_node_alt_state_exp_hphobes.size(); ++jj ) {
				//TR << ii_node_alt_state_exp_hphobes[ jj ] << ", ";
				TS_ASSERT( ( rotamer->atom_type( ii_node_alt_state_exp_hphobes[ jj ] ).element() == "C" ) ||
					( rotamer->atom_type( ii_node_alt_state_exp_hphobes[ jj ] ).element() == "S" ) );
				TS_ASSERT( pdhig->get_hpatch_node( ii )->get_alternate_state_sasa( ii_node_alt_state_exp_hphobes[ jj ] ) > 0.0 );
			}
			//TR << "]" << std::endl;
		}

		// BG nodes, current state and alt state
		for ( Size ii = 1; ii <= (Size)pdhig->get_num_background_nodes(); ++ii ) {
			utility::vector1< Size > const & ii_bg_node_curr_state_exp_hphobes = pdhig->get_hpatch_bg_node( ii )->curr_state_exp_hphobes();
			conformation::ResidueCOP rotamer = pdhig->get_hpatch_bg_node( ii )->get_rotamer();
			//TR << "test_curr_state_exp_hphobes(): node " << ii << ": ii_bg_node_curr_state_exp_hphobes: [ ";
			for ( Size jj=1; jj <= ii_bg_node_curr_state_exp_hphobes.size(); ++jj ) {
				//TR << ii_bg_node_curr_state_exp_hphobes[ jj ] << ", ";
				TS_ASSERT( ( rotamer->atom_type( ii_bg_node_curr_state_exp_hphobes[ jj ] ).element() == "C" ) ||
					( rotamer->atom_type( ii_bg_node_curr_state_exp_hphobes[ jj ] ).element() == "S" ) );
				TS_ASSERT( pdhig->get_hpatch_bg_node( ii )->get_current_sasa( ii_bg_node_curr_state_exp_hphobes[ jj ] ) > 0.0 );
			}
			//TR << "], ";

			utility::vector1< Size > const & ii_bg_node_alt_state_exp_hphobes = pdhig->get_hpatch_bg_node( ii )->alt_state_exp_hphobes();
			//TR << "test_curr_state_exp_hphobes(): ii_bg_node_alt_state_exp_hphobes: [ ";
			for ( Size jj=1; jj <= ii_bg_node_alt_state_exp_hphobes.size(); ++jj ) {
				//TR << ii_bg_node_alt_state_exp_hphobes[ jj ] << ", ";
				TS_ASSERT( ( rotamer->atom_type( ii_bg_node_alt_state_exp_hphobes[ jj ] ).element() == "C" ) ||
					( rotamer->atom_type( ii_bg_node_alt_state_exp_hphobes[ jj ] ).element() == "S" ) );
				TS_ASSERT( pdhig->get_hpatch_bg_node( ii )->get_alternate_sasa( ii_bg_node_alt_state_exp_hphobes[ jj ] ) > 0.0 );
			}
			//TR << "]" << std::endl;
		}

	}


	/// @details
	/// Tests the function consider_substitution() (hereafter cs()).
	/// cs() takes a position and a new state and updates the SASA for that node and all its neighboring nodes/bgnodes.
	/// It's possible that cs() will procrastinate the calculation if the change in the PD energy is really bad, but
	/// the IG shouldn't procrastinate anything after just being prepped for simA and init'd.
	///
	void test_consider_substitution() {
		TR << "Running test_consider_substitution..." << std::endl;

		// In this new version of the IG, every substitution will lead to a change in the score.

		// P->P: commit state 382(ARG) on MR 5 (PDB: 8, LYS)
		// this sub only produces SASA changes on bg node indexes 4 and 5 (1-3 don't change)
		pdhig->consider_substitution( 5, 382, delta_energy, previous_energy_for_node );
		TS_ASSERT_DELTA( delta_energy, 0.7064, 0.05 );
		pdhig->commit_considered_substitution();

		// HP->P: commit state 186 (TYR) on MR 2 (GLN)
		pdhig->consider_substitution( 2, 186, delta_energy, previous_energy_for_node );
		TS_ASSERT_DELTA( delta_energy, 0.7256, 0.05 );
		pdhig->commit_considered_substitution();

		int states_a[] = { 0, 186, 0, 0, 382, 0 };
		const size_t size = sizeof(states_a) / sizeof (states_a[0]);
		std::vector<int> correct_state_a(states_a, states_a + size);
		TS_ASSERT_EQUALS( pdhig->get_network_state(), correct_state_a );

		/*for ( Size ii=1; ii <= (Size)pdhig->get_num_nodes(); ++ii ) {
		TR << "node " << I(2,ii) << ": " << std::endl;
		((pdhig->get_hpatch_node(ii))->get_current_state_rotamer_dots()).print( std::cout );
		}
		for ( Size ii=1; ii <= (Size)pdhig->get_num_background_nodes(); ++ii ) {
		TR << "bgnode " << I(2,ii) << ": " << std::endl;
		((pdhig->get_hpatch_bg_node(ii))->get_current_state_rotamer_dots()).print( std::cout );
		}*/

		// need to also test here that the nodes have updated their hASAs correctly!
		utility::vector1< Real > node_sasas( (Size)pdhig->get_num_nodes(), 0.0 );
		utility::vector1< Real > bgnode_sasas ( (Size)pdhig->get_num_background_nodes(), 0.0 );
		pdhig->get_all_sasas( node_sasas, bgnode_sasas );

		Real correct_node_sasas[] = { 0.0, 186.0778, 0.0, 0.0, 328.7393, 0.0 };
		Real correct_bgnode_sasas[] = { 233.5590, 108.5452, 213.1477, 63.8602, 66.3769, 140.9998, 151.6156, 19.0419, 101.2664, 151.2777, 118.2961, 73.8983, 37.8562, 240.0214 };


		// all the Nodes except 2 and 5 should have 0.0 for their SASA because they should still be in the unassigned state
		for ( Size ii=1; ii <= node_sasas.size(); ++ii ) {
			TS_ASSERT_DELTA( node_sasas[ ii ], correct_node_sasas[ ii-1 ], 0.01 );
		}
		// all the background nodes should have some value for SASA
		for ( Size ii=1; ii <= bgnode_sasas.size(); ++ii ) {
			TS_ASSERT_DELTA( bgnode_sasas[ ii ], correct_bgnode_sasas[ ii-1 ], 0.01 );
		}

		// HP(a)->HP(b): commit state 61 PHE on TRP-6 (MR: 3) #ha: 160
		pdhig->consider_substitution( 3, 61, delta_energy, previous_energy_for_node );
		TS_ASSERT_DELTA( delta_energy, -0.8923, 0.05 );
		currentenergy = pdhig->commit_considered_substitution();

		// HP(a)->P(b): commit state 25 GLU on LEU-7 (MR: 4) #ha: 161
		pdhig->consider_substitution( 4, 25, delta_energy, previous_energy_for_node );
		TS_ASSERT_DELTA( delta_energy, -1.2401, 0.05 );
		currentenergy = pdhig->commit_considered_substitution();

		// HP->P: commit state 435 ARG on ILE-4 (MR: 1) #ha: 164
		pdhig->consider_substitution( 1, 435, delta_energy, previous_energy_for_node );
		TS_ASSERT_DELTA( delta_energy, -2.01762, 0.05 );
		currentenergy = pdhig->commit_considered_substitution();

		// Check for correct SASAs on all nodes again
		pdhig->get_all_sasas( node_sasas, bgnode_sasas );

		Real correct_node_sasas_b[] = { 191.2383, 162.4853, 6.1956, 116.6548, 250.0062, 0.0 };
		Real correct_bgnode_sasas_b[] = { 204.0776, 104.0648, 135.6771, 60.4861, 39.9076, 114.8925, 151.6156, 19.0419, 101.2664, 141.8631, 113.8156, 72.6259, 29.9482, 240.0214 };

		for ( Size ii=1; ii <= node_sasas.size(); ++ii ) {
			TS_ASSERT_DELTA( node_sasas[ ii ], correct_node_sasas_b[ ii-1 ], 0.01 );
		}
		// all the background nodes should have some value for SASA
		for ( Size ii=1; ii <= bgnode_sasas.size(); ++ii ) {
			TS_ASSERT_DELTA( bgnode_sasas[ ii ], correct_bgnode_sasas_b[ ii-1 ], 0.01 );
		}

	}


	/// @details
	/// Tests the function calculate_alt_state_hpatch_score(). Have to make 6 commits minimum before the score gets evaluated
	/// because the score only gets calculated if all the nodes are in some state (not the unassigned state).
	///
	/// This test will also test the methods register_fc_node_affected_by_rotsub and register_bg_node_affected_by_rotsub
	/// as these functions also only really do something when there are no nodes left unassigned.
	///
	void test_calculate_alt_state_hpatch_score() {
		TR << "Running test_calculate_alt_state_hpatch_score..." << std::endl;

		// In this new version of the SIG, every substitution will lead to a change in the score.

		// P->P: commit state 382(ARG) on MR 5 (PDB: 8, LYS)
		// this sub only produces SASA changes on bg node indexes 4 and 5 (1-3 don't change)
		pdhig->consider_substitution( 5, 382, delta_energy, previous_energy_for_node );
		pdhig->commit_considered_substitution();

		// HP->P: commit state 186 (TYR) on MR 2 (GLN)
		pdhig->consider_substitution( 2, 186, delta_energy, previous_energy_for_node );
		pdhig->commit_considered_substitution();

		// HP(a)->HP(b): commit state 61 PHE on TRP-6 (MR: 3) #ha: 160
		pdhig->consider_substitution( 3, 61, delta_energy, previous_energy_for_node );
		pdhig->commit_considered_substitution();

		// HP(a)->P(b): commit state 25 GLU on LEU-7 (MR: 4) #ha: 161
		pdhig->consider_substitution( 4, 25, delta_energy, previous_energy_for_node );
		pdhig->commit_considered_substitution();

		// HP->P: commit state 435 ARG on ILE-4 (MR: 1) #ha: 164
		pdhig->consider_substitution( 1, 435, delta_energy, previous_energy_for_node );
		pdhig->commit_considered_substitution();

		// commit state 17 SER on GLY-11 (MR: 6)
		pdhig->consider_substitution( 6, 17, delta_energy, previous_energy_for_node );
		pdhig->commit_considered_substitution();

		// all the Nodes are in an assigned state now. considering one more sub should cause a score evaluation.

		// commit state 551(TRP) on MR 5 (PDB: 8, LYS)
		pdhig->consider_substitution( 5, 551, delta_energy, previous_energy_for_node );
		//currentenergy = pdhig->commit_considered_substitution();

		// make sure energy is correct
		TS_ASSERT_DELTA( delta_energy, 4.1270, 0.01 );

		int states_a[] = { 435, 186, 61, 25, 382, 17 };
		const size_t size = sizeof(states_a) / sizeof (states_a[0]);
		std::vector<int> correct_state_a(states_a, states_a + size);
		TS_ASSERT_EQUALS( pdhig->get_network_state(), correct_state_a );
	}


	/// @details
	/// The main things to test with commit sub() is that the total energy returned is correct and that the node counts (for
	/// the changing node *and* all neighboring nodes) are updated.  It's possible that consider doesn't actually compute
	/// the correct energy because of computation procrastination.  That would also mean the node counts would be inaccurate
	/// until after the commit occurred.
	///
	void test_commit_substitution() {
		TR << "Running test_commit_substitution..." << std::endl;

		// Tests we'll have here:
		// 1) whether the delta_energy returned by the method is correct for the case where the calculation is not
		// procrastinated, and 2) when it is procrastinated.
		// Another thing to test for is that the hASAs reset if one sub is not commit'd but a following sub is.

		// P->HP: commit state 186 (TYR) on MR 2 (GLN)
		pdhig->consider_substitution( 2, 186, delta_energy, previous_energy_for_node );
		TS_ASSERT_DELTA( delta_energy, 1.12343, 0.05 );
		pdhig->commit_considered_substitution();

		int states_a[] = { 0, 186, 0, 0, 0, 0 };
		const size_t size_a = sizeof(states_a) / sizeof (states_a[0]);
		std::vector<int> correct_state_a(states_a, states_a + size_a);
		TS_ASSERT_EQUALS( pdhig->get_network_state(), correct_state_a );

		// need to also test here that the nodes have updated their hASAs correctly!
		utility::vector1< Real > node_sasas( (Size)pdhig->get_num_nodes(), 0.0 );
		utility::vector1< Real > bgnode_sasas ( (Size)pdhig->get_num_background_nodes(), 0.0 );
		pdhig->get_all_sasas( node_sasas, bgnode_sasas );

		Real correct_node_sasas[] = { 0.0, 213.5716, 0.0, 0.0, 0.0, 0.0 };
		Real correct_bgnode_sasas[] = { 233.5590, 108.5452, 227.7092, 110.3281, 82.9174, 140.9998, 151.6156, 19.0419, 101.2664, 151.2777, 118.2961, 73.8983, 37.8562, 240.0214 };

		// all the Nodes except 2 and 5 should have 0.0 for their SASA because they should still be in the unassigned state
		for ( Size ii=1; ii <= node_sasas.size(); ++ii ) {
			TS_ASSERT_DELTA( node_sasas[ ii ], correct_node_sasas[ ii-1 ], 0.01 );
		}
		// all the background nodes should have some value for SASA
		for ( Size ii=1; ii <= bgnode_sasas.size(); ++ii ) {
			TS_ASSERT_DELTA( bgnode_sasas[ ii ], correct_bgnode_sasas[ ii-1 ], 0.01 );
		}

		// HP->P: commit state 435 ARG on ILE-4 (MR: 1)
		pdhig->consider_substitution( 1, 435, delta_energy, previous_energy_for_node );
		// don't commit

		// HP(a)->HP(b): commit state 25 PHE on LEU-7 (MR: 4)
		pdhig->consider_substitution( 4, 25, delta_energy, previous_energy_for_node );
		TS_ASSERT_DELTA( delta_energy, -0.735305, 0.05 );
		currentenergy = pdhig->commit_considered_substitution();

		int states_b[] = { 0, 186, 0, 25, 0, 0 };
		const size_t size_b = sizeof(states_b) / sizeof (states_b[0]);
		std::vector<int> correct_state_b(states_b, states_b + size_b);
		TS_ASSERT_EQUALS( pdhig->get_network_state(), correct_state_b );


		// --- make a really bad mutation which should procrastinate the surface calculation
		pdhig->set_observed_sufficient_boolean_true(); // this function call ensures that the surface calculation is procrastinated

		pdhig->consider_substitution( 6, 16, delta_energy, previous_energy_for_node );
		// delta_energy should be an estimated value (equal to the PD deltaE alone) if the graph is procrastinating correctly
		TS_ASSERT_DELTA( delta_energy, 192.493, 5.0 );

		// to make sure the IG procrastinated, check the node state before and after the commit
		int states_c[] = { 0, 186, 0, 25, 0, 0 };
		const size_t size_c = sizeof(states_c) / sizeof (states_c[0]);
		std::vector<int> correct_state_c(states_c, states_c + size_c);
		TS_ASSERT_EQUALS( pdhig->get_network_state(), correct_state_c );

		currentenergy = pdhig->commit_considered_substitution();
		TS_ASSERT_DELTA( currentenergy, 192.8810, 5.0 );

		int states_d[] = { 0, 186, 0, 25, 0, 16 };
		std::vector<int> correct_state_d(states_d, states_d + size_c);
		TS_ASSERT_EQUALS( pdhig->get_network_state(), correct_state_d );

	}

	/// @details
	/// Near the end of sims, lots of rotamers are tried (which change the alt state counts) but then aren't committed.
	/// This test ensures that the graph is resetting state correctly in those cases.
	/// blanket_reset is implicitly tested in the unit tests above...
	///
	//void x_test_blanket_reset_alt_state_counts() {
	// TR << "Running test_blanket_reset_alt_state_counts..." << std::endl;
	//}

	/// @details
	/// Make some random commits and make sure the total energy current state assignment method is returning the same
	/// thing that was computed in commit_sub().
	///
	void test_get_energy_current_state_assignment() {
		TR << "Running test_get_energy_current_state_assignment..." << std::endl;

		pdhig->consider_substitution( 2, 106, delta_energy, previous_energy_for_node );
		currentenergy = pdhig->commit_considered_substitution();
		TS_ASSERT( currentenergy == pdhig->get_energy_current_state_assignment() );

		pdhig->consider_substitution( 5, 10, delta_energy, previous_energy_for_node );
		currentenergy = pdhig->commit_considered_substitution();
		TS_ASSERT( currentenergy == pdhig->get_energy_current_state_assignment() );

		pdhig->consider_substitution( 6, 20, delta_energy, previous_energy_for_node );
		currentenergy = pdhig->commit_considered_substitution();
		TS_ASSERT( currentenergy == pdhig->get_energy_current_state_assignment() );

	}


	/// @brief
	/// a simple packing run that uses surface scoring and a standard PD interaction graph
	///
	void test_partial_redesign_using_pd_ig() {
		TR << "Running test_partial_redesign_using_pd_ig..." << std::endl;

		// Set these really cool options Mike Tyka added to reduce the number of cycles the annealer runs
		// so that this test doesn't take forever.
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		basic::options::option[ basic::options::OptionKeys::packing::outeriterations_scaling ].value( 0.2 );
		basic::options::option[ basic::options::OptionKeys::packing::inneriterations_scaling ].value( 0.05 );

		// --- InteractionGraph ---
		bool start_with_current = false;
		bool calc_rot_freq = false;
		ObjexxFCL::FArray1D_int current_rot_index; current_rot_index.dimension( pose.total_residue(), 0 );
		ObjexxFCL::FArray1D< PackerEnergy > rot_freq; rot_freq.dimension( pdhig->get_num_total_states(), 0.0 );
		utility::vector0<int> rot_to_pack;

		annealer::SimAnnealerBaseOP redesign_annealer = annealer::AnnealerFactory::create_annealer(
			designtask, rot_to_pack, bestrotamer_at_seqpos, bestenergy, start_with_current, utility::pointer::dynamic_pointer_cast<interaction_graph::AnnealableGraphBase>(pdhig), rotsets, current_rot_index, calc_rot_freq, rot_freq );

		// temperature isn't so important, but to make things easy use the SA setup_temp() method
		ObjexxFCL::FArray1D_float loopenergy( 500, 0.0 );  // hardcore the number of loops for this array to the maxnumberofouteriterations
		annealer->setup_temperature( loopenergy, 1 );  // 1 would be the first iteration of outer loop

		redesign_annealer->run();

		// check to make sure the SASAs are all correct using sasa.cc
		// after the annealer is done runnning, place new rotamers on the input pose
		for ( core::uint ii = 1; ii <= rotsets->nmoltenres(); ++ii ) {
			core::uint iiresid = rotsets->moltenres_2_resid( ii );
			core::uint iibestrot = rotsets->rotid_on_moltenresidue( bestrotamer_at_seqpos( iiresid ) );
			conformation::ResidueCOP bestrot( rotsets->rotamer_set_for_moltenresidue( ii )->rotamer( iibestrot ) );

			conformation::ResidueOP newresidue( bestrot->create_residue() );
			pose.replace_residue( iiresid, *newresidue, false );
		}

		id::AtomID_Map< Real > atom_sasa;
		utility::vector1< Real > rsd_sasa;
		Real probe_radius = 1.4;

		// create an atom_subset mask that will tell calc_per_atom_sasa to ignore all H's
		id::AtomID_Map< bool > atom_subset;
		atom_subset.resize( pose.n_residue() );
		// init all heavy atoms to true and all H's to false
		for ( Size ii=1; ii <= pose.n_residue(); ++ii ) {
			atom_subset.resize( ii, pose.residue_type(ii).natoms(), false );
			for ( Size jj = 1; jj <= pose.residue_type(ii).nheavyatoms(); ++jj ) {
				atom_subset[ ii ][ jj ] = true;
			}
		}

		core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius, false /* no big polar H */,
			atom_subset, true /* use_naccess_sasa_radii */, true /* expand polar radii */ );

		// put all the Nodes/BGNodes into the final best state, not the state that simA ended at
		ObjexxFCL::FArray1D_int states;
		states.dimension( rotsets->nmoltenres() );
		for ( core::uint ii = 1; ii <= rotsets->nmoltenres(); ++ii ) {
			core::uint iiresid = rotsets->moltenres_2_resid( ii );
			core::uint iibestrot = rotsets->rotid_on_moltenresidue( bestrotamer_at_seqpos( iiresid ) );
			states( ii ) = iibestrot;
		}
		pdhig->set_network_state( states );

		utility::vector1< Real > node_sasas( (Size)pdhig->get_num_nodes(), 0.0 );
		utility::vector1< Real > bgnode_sasas ( (Size)pdhig->get_num_background_nodes(), 0.0 );
		pdhig->get_all_sasas( node_sasas, bgnode_sasas );

		for ( Size ii=1; ii <= node_sasas.size(); ++ii ) {
			TS_ASSERT_DELTA( node_sasas[ ii ], rsd_sasa[ rotsets->moltenres_2_resid( ii ) ], 0.01 );
		}
		// all the background nodes should have some value for SASA
		for ( Size ii=1; ii <= bgnode_sasas.size(); ++ii ) {
			TS_ASSERT_DELTA( bgnode_sasas[ ii ], rsd_sasa[ pdhig->bg_node_2_resid( ii ) ], 0.01 );
		}

	}

	/// @brief
	/// A unit (maybe more of an integration test) that does a short design run using the surface score and a linear memory
	/// interaction graph.
	///
	void test_partial_redesign_using_linmem_ig() {
		TR << "Running test_partial_redesign_using_linmem_ig..." << std::endl;

		// Set these really cool options Mike Tyka added to reduce the number of cycles the annealer runs
		// so that this test doesn't take forever.
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		basic::options::option[ basic::options::OptionKeys::packing::outeriterations_scaling ].value( 0.2 );
		basic::options::option[ basic::options::OptionKeys::packing::inneriterations_scaling ].value( 0.05 );

		// setting an option in the code, very bad coding practice. there's no way to set the value of the recent
		// history size via a PackerTask. only way that exists now is to use the command line which we can't do for
		// just this one test and not the other tests in this Suite. so manually set the option here. it's just a unit test!
		//basic::options::option[ basic::options::OptionKeys::packing::linmem_ig ].value( 10 );

		// for this test, we don't want to use the PDIG that's created in the fixture

		// TR << "Instantiating LinearMemoryHPatchInteractionGraph..." << std::endl;
		lmhig = interaction_graph::LinearMemoryHPatchInteractionGraphOP( new interaction_graph::LinearMemoryHPatchInteractionGraph( designtask->num_to_be_packed() ) );
		lmhig->set_pose( pose );
		lmhig->set_packer_task( *designtask );
		lmhig->set_score_function( *scorefxn );
		lmhig->set_rotamer_sets( *rotsets );

		// lmhig is a LinearMem HIG OP; compute_energies() wants an IGBase OP. Need to dereference the
		// OP to get the HIG and then cast to an IGBase.

		TR << "\tcomputing rotamer pair energies..." << std::endl;
		// compute_energies() does some initialization of the interaction graph and computes the energies
		rotsets->compute_energies( pose, *scorefxn, packer_neighbor_graph, static_cast< interaction_graph::InteractionGraphBaseOP >(lmhig) );

		/// Parameters passed by reference in annealers constructor to which it writes at the completion of sim annealing.
		bool start_with_current = false;
		bool calc_rot_freq = false;
		ObjexxFCL::FArray1D_int current_rot_index; current_rot_index.dimension( pose.total_residue(), 0 );
		ObjexxFCL::FArray1D< PackerEnergy > linmem_ig_test_rot_freq( lmhig->get_num_total_states(), 0.0 );
		utility::vector0<int> rot_to_pack;

		annealer::SimAnnealerBaseOP redesign_annealer = annealer::AnnealerFactory::create_annealer(
			designtask, rot_to_pack, bestrotamer_at_seqpos, bestenergy, start_with_current, utility::pointer::dynamic_pointer_cast<interaction_graph::AnnealableGraphBase>(lmhig),
			rotsets, current_rot_index, calc_rot_freq, linmem_ig_test_rot_freq );

		redesign_annealer->run();

		// check to make sure the SASAs are all correct using sasa.cc
		// after the annealer is done runnning, place new rotamers on the input pose
		for ( core::uint ii = 1; ii <= rotsets->nmoltenres(); ++ii ) {
			core::uint iiresid = rotsets->moltenres_2_resid( ii );
			core::uint iibestrot = rotsets->rotid_on_moltenresidue( bestrotamer_at_seqpos( iiresid ) );
			conformation::ResidueCOP bestrot( rotsets->rotamer_set_for_moltenresidue( ii )->rotamer( iibestrot ) );

			conformation::ResidueOP newresidue( bestrot->create_residue() );
			pose.replace_residue( iiresid, *newresidue, false );
		}

		id::AtomID_Map< Real > atom_sasa;
		utility::vector1< Real > rsd_sasa;
		Real probe_radius = 1.4;

		// create an atom_subset mask that will tell calc_per_atom_sasa to ignore all H's
		id::AtomID_Map< bool > atom_subset;
		atom_subset.resize( pose.n_residue() );
		// init all heavy atoms to true and all H's to false
		for ( Size ii=1; ii <= pose.n_residue(); ++ii ) {
			atom_subset.resize( ii, pose.residue_type(ii).natoms(), false );
			for ( Size jj = 1; jj <= pose.residue_type(ii).nheavyatoms(); ++jj ) {
				atom_subset[ ii ][ jj ] = true;
			}
		}

		core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius, false /* no big polar H */,
			atom_subset, true /* use_naccess_sasa_radii */, true /* expand polar radii */ );

		// put all the Nodes/BGNodes into the final best state, not the state that simA ended at
		ObjexxFCL::FArray1D_int states;
		states.dimension( rotsets->nmoltenres() );
		for ( core::uint ii = 1; ii <= rotsets->nmoltenres(); ++ii ) {
			core::uint iiresid = rotsets->moltenres_2_resid( ii );
			core::uint iibestrot = rotsets->rotid_on_moltenresidue( bestrotamer_at_seqpos( iiresid ) );
			states( ii ) = iibestrot;
		}
		lmhig->set_network_state( states );

		utility::vector1< Real > node_sasas( (Size)lmhig->get_num_nodes(), 0.0 );
		utility::vector1< Real > bgnode_sasas ( (Size)lmhig->get_num_background_nodes(), 0.0 );
		lmhig->get_all_sasas( node_sasas, bgnode_sasas );

		for ( Size ii=1; ii <= node_sasas.size(); ++ii ) {
			TS_ASSERT_DELTA( node_sasas[ ii ], rsd_sasa[ rotsets->moltenres_2_resid( ii ) ], 0.01 );
		}
		// all the background nodes should have some value for SASA
		for ( Size ii=1; ii <= bgnode_sasas.size(); ++ii ) {
			TS_ASSERT_DELTA( bgnode_sasas[ ii ], rsd_sasa[ lmhig->bg_node_2_resid( ii ) ], 0.01 );
		}

	}

};
