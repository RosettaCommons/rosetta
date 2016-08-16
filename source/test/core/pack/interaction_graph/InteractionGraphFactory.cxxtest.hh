// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/InteractionGraphFactory.cxxtest.hh
/// @brief  test for the InteractionGraphFactory
/// @author Steven Lewis

// Test framework headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Core Headers
#include <core/pack/interaction_graph/DensePDInteractionGraph.hh>
#include <core/pack/interaction_graph/DoubleLazyInteractionGraph.hh>
#include <core/pack/interaction_graph/InteractionGraphFactory.hh>
#include <core/pack/interaction_graph/LazyInteractionGraph.hh>
#include <core/pack/interaction_graph/LinearMemoryInteractionGraph.hh>
#include <core/pack/interaction_graph/PDInteractionGraph.hh>
#include <core/pack/interaction_graph/SurfaceInteractionGraph.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/pack/packer_neighbors.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>

#include <core/util/SwitchResidueTypeSet.hh>

//Auto Headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1.io.hh>


//using namespace core;
using namespace core::pack;
using namespace core::pack::interaction_graph;

// --------------- Test Class --------------- //
/// @details This suite of tests covers the InteractionGraphFactory.  That class's job is to examine the option system/PackerTask and determine the appropriate flavor of InteractionGraph for the packing at hand.  This test creates environments, runs the factory, and uses dynamic_cast to ensure that the created IG is of the correct type.
class InteractionGraphFactoryTests : public CxxTest::TestSuite {

public:

	bool suite_initialized;

	// Shared data elements go here.
	core::pose::Pose pose;
	rotamer_set::RotamerSetsOP rotsets;
	core::scoring::ScoreFunctionOP scorefxn;
	core::graph::GraphOP packer_neighbor_graph;
	task::PackerTaskOP packertask;


	// --------------- Suite-level Fixture --------------- //

	InteractionGraphFactoryTests() {

		core_init();
		pose = create_twores_1ubq_pose();

		// --- ScoreFunction ---
		// create a score function using the standard packer weights
		scorefxn = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction() );
		//scorefxn->set_weight( scoring::surface, 0.5 );


	}

	virtual ~InteractionGraphFactoryTests() {}

	static InteractionGraphFactoryTests *createSuite() {
		return new InteractionGraphFactoryTests();
	}

	static void destroySuite( InteractionGraphFactoryTests *suite ) {
		delete suite;
	}

	// --------------- Test Fixture --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case.

	void setUp() {
	}

	// Shared finalization goes here.
	void tearDown() {}

	//utility function - we don't care about any of these rotamers but we have to do it anyway for the IGFactory
	//relies on class member variables instead of passing stuff around
	void make_rotset(){
		rotsets = rotamer_set::RotamerSetsOP( new rotamer_set::RotamerSets() );
		rotsets->set_task( packertask );
		packer_neighbor_graph = create_packer_graph( pose, *scorefxn, packertask );
		rotsets->build_rotamers( pose, *scorefxn, packer_neighbor_graph );
		rotsets->prepare_sets_for_packing( pose, *scorefxn );
		return;
	}


public:

	// --------------- Test Cases --------------- //

	//As of this writing, we have the following options:

	//If linmem_if is set in the PackerTask
	//LinearMemoryInteractionGraph
	//LinearMemorySurfaceInteractionGraph
	//    if surface weight !=0 and and we are designing

	//else if we are designing, altering >=1 residue, it has >=1 rotamers, and it's not centroid
	//PDInteractionGraph
	//PDSurfaceInteractionGraph
	//    if surface weight !=0, and we are designing
	//LazyInteractionGraph
	//    if lazy_ig is set in the PackerTask
	//DoubleLazyInteractionGraph
	//    if double_lazy_if is set in the PackerTask

	//else (meaning, no linmem_ig, not designing, in centroid mode, or we have no rotamers)
	//DensePDInteractionGraph

	/// @details we should get a linmemIG if linmem_ig is set in the packer task and surface weight is 0 in the scorefunction
	void test_LinearMemoryInteractionGraph() {
		//surface weight 0
		scorefxn->set_weight( core::scoring::surface, 0 );
		(*scorefxn)(pose);

		//create packer task, set linmem_ig
		packertask = core::pack::task::TaskFactory::create_packer_task( pose );
		packertask->or_linmem_ig(true);

		//create the rotamer sets
		make_rotset();

		//TEST that the IGFactory returns a linmem_ig
		TS_ASSERT( dynamic_cast<LinearMemoryInteractionGraph*>( (InteractionGraphFactory::create_interaction_graph( *packertask, *rotsets, pose, *scorefxn ).get()) ));
		//let's test a random fail state for kicks
		TS_ASSERT( !(dynamic_cast<PDInteractionGraph*>( (InteractionGraphFactory::create_interaction_graph( *packertask, *rotsets, pose, *scorefxn ).get()) )) );

		return;
	}

	/// @brief we should get a linmemsurfaceIG if linmem_ig is set in the packer task and surface weight is not 0 in the scorefunction
	void test_LinearMemorySurfaceInteractionGraph() {
		//surface weight 1
		scorefxn->set_weight( core::scoring::surface, 1 );
		(*scorefxn)(pose);

		//create packer task, set linmem_ig
		packertask = core::pack::task::TaskFactory::create_packer_task( pose );
		packertask->or_linmem_ig(true);

		//create the rotamer sets
		make_rotset();

		//TEST what the IGFactory returns
		TS_ASSERT( dynamic_cast<LinearMemorySurfaceInteractionGraph*>( (InteractionGraphFactory::create_interaction_graph( *packertask, *rotsets, pose, *scorefxn ).get()) ));
		//let's test a random fail state for kicks
		TS_ASSERT( !(dynamic_cast<PDInteractionGraph*>( (InteractionGraphFactory::create_interaction_graph( *packertask, *rotsets, pose, *scorefxn ).get()) )) );

		return;
	}

	/// @brief we should get this if we are designing, have rotamers, aren't centroid, and surface is nonzero
	void test_PDSurfaceInteractionGraph() {
		//surface weight 1
		scorefxn->set_weight( core::scoring::surface, 1 );
		(*scorefxn)(pose);

		//create packer task
		packertask = core::pack::task::TaskFactory::create_packer_task( pose );

		//create the rotamer sets
		make_rotset();

		//TEST what the IGFactory returns
		TS_ASSERT( dynamic_cast<PDSurfaceInteractionGraph*>( (InteractionGraphFactory::create_interaction_graph( *packertask, *rotsets, pose, *scorefxn ).get()) ));
		//let's test a random fail state for kicks
		TS_ASSERT( !(dynamic_cast<LinearMemoryInteractionGraph*>( (InteractionGraphFactory::create_interaction_graph( *packertask, *rotsets, pose, *scorefxn ).get()) )) );

		return;
	}

	/// @brief we should get this if we are designing, have rotamers, aren't centroid, and surface is zero
	void test_PDInteractionGraph() {
		//surface weight 0
		scorefxn->set_weight( core::scoring::surface, 0 );
		(*scorefxn)(pose);

		//create packer task
		packertask = core::pack::task::TaskFactory::create_packer_task( pose );

		//create the rotamer sets
		make_rotset();

		//TEST what the IGFactory returns
		TS_ASSERT( dynamic_cast<PDInteractionGraph*>( (InteractionGraphFactory::create_interaction_graph( *packertask, *rotsets, pose, *scorefxn ).get()) ));
		//let's test a random fail state for kicks
		TS_ASSERT( !(dynamic_cast<LinearMemoryInteractionGraph*>( (InteractionGraphFactory::create_interaction_graph( *packertask, *rotsets, pose, *scorefxn ).get()) )) );

		return;
	}

	/// @brief we should get this if we are designing, have rotamers, aren't centroid, and surface is zero, and lazy_ig is set in the packertask
	void test_LazyInteractionGraph() {
		//surface weight 0
		scorefxn->set_weight( core::scoring::surface, 0 );
		(*scorefxn)(pose);

		//create packer task
		packertask = core::pack::task::TaskFactory::create_packer_task( pose );
		packertask->or_lazy_ig(true);

		//create the rotamer sets
		make_rotset();

		//TEST what the IGFactory returns
		TS_ASSERT( dynamic_cast<LazyInteractionGraph*>( (InteractionGraphFactory::create_interaction_graph( *packertask, *rotsets, pose, *scorefxn ).get()) ));
		//let's test a random fail state for kicks
		TS_ASSERT( !(dynamic_cast<LinearMemoryInteractionGraph*>( (InteractionGraphFactory::create_interaction_graph( *packertask, *rotsets, pose, *scorefxn ).get()) )) );

		return;
	}

	/// @brief we should get this if we are designing, have rotamers, aren't centroid, and surface is zero, and double_lazy_ig is set in the packertask
	void test_DoubleLazyInteractionGraph() {
		//surface weight 0
		scorefxn->set_weight( core::scoring::surface, 0 );
		(*scorefxn)(pose);

		//create packer task
		packertask = core::pack::task::TaskFactory::create_packer_task( pose );
		packertask->or_double_lazy_ig(true);

		//create the rotamer sets
		make_rotset();

		//TEST what the IGFactory returns
		TS_ASSERT( dynamic_cast<DoubleLazyInteractionGraph*>( (InteractionGraphFactory::create_interaction_graph( *packertask, *rotsets, pose, *scorefxn ).get()) ));
		//let's test a random fail state for kicks
		TS_ASSERT( !(dynamic_cast<LinearMemoryInteractionGraph*>( (InteractionGraphFactory::create_interaction_graph( *packertask, *rotsets, pose, *scorefxn ).get()) )) );

		return;
	}

	/// @brief we should get this if we aren't designing or are centroid, and linmem_ig is not set in the packertask
	void test_DensePDInteractionGraph() {
		//TEST 1: not designing
		//surface weight 0
		scorefxn->set_weight( core::scoring::surface, 0 );
		(*scorefxn)(pose);

		//create packer task
		packertask = core::pack::task::TaskFactory::create_packer_task( pose );
		packertask->restrict_to_repacking();

		//create the rotamer sets
		make_rotset();

		//TEST what the IGFactory returns
		TS_ASSERT( dynamic_cast<DensePDInteractionGraph*>( (InteractionGraphFactory::create_interaction_graph( *packertask, *rotsets, pose, *scorefxn ).get()) ));
		//let's test a random fail state for kicks
		TS_ASSERT( !(dynamic_cast<LinearMemoryInteractionGraph*>( (InteractionGraphFactory::create_interaction_graph( *packertask, *rotsets, pose, *scorefxn ).get()) )) );


		//TEST 2: no rotamers
		//surface weight 0
		scorefxn->set_weight( core::scoring::surface, 0 );
		(*scorefxn)(pose);

		//create packer task
		packertask = core::pack::task::TaskFactory::create_packer_task( pose );
		for ( core::Size i(1), end(packertask->total_residue()); i<=end; ++i ) {
			packertask->nonconst_residue_task(i).prevent_repacking();
		}

		//create the rotamer sets
		make_rotset();

		//TEST what the IGFactory returns
		TS_ASSERT( dynamic_cast<DensePDInteractionGraph*>( (InteractionGraphFactory::create_interaction_graph( *packertask, *rotsets, pose, *scorefxn ).get()) ));
		//let's test a random fail state for kicks
		TS_ASSERT( !(dynamic_cast<LinearMemoryInteractionGraph*>( (InteractionGraphFactory::create_interaction_graph( *packertask, *rotsets, pose, *scorefxn ).get()) )) );

		//TEST 3: CENTROID!  whoo
		core::util::switch_to_residue_type_set(pose, core::chemical::CENTROID);

		//surface weight 0
		scorefxn->set_weight( core::scoring::surface, 0 );
		(*scorefxn)(pose);

		//create packer task
		packertask = core::pack::task::TaskFactory::create_packer_task( pose );

		//create the rotamer sets
		make_rotset();

		//TEST what the IGFactory returns
		TS_ASSERT( dynamic_cast<DensePDInteractionGraph*>( (InteractionGraphFactory::create_interaction_graph( *packertask, *rotsets, pose, *scorefxn ).get()) ));
		//let's test a random fail state for kicks
		TS_ASSERT( !(dynamic_cast<LinearMemoryInteractionGraph*>( (InteractionGraphFactory::create_interaction_graph( *packertask, *rotsets, pose, *scorefxn ).get()) )) );

		return;
	}

};

