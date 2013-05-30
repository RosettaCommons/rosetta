// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/FASTERInteractionGraph.cxxtest.hh
/// @brief  test suite for the FASTER interaction graph
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test framework headers
#include <cxxtest/TestSuite.h>

// Core Headers
#include <core/pack/interaction_graph/FASTERInteractionGraph.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>

#include <core/graph/Graph.hh>

// AUTO-REMOVED #include <core/pose/Pose.hh>
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


class FASTERInteractionGraphTests : public CxxTest::TestSuite {
public:

	void setUp() {
		core_init_with_additional_options( "-restore_pre_talaris_2013_behavior" );
	}

	void test_instantiate_FASTER_ig() {
		using namespace core::chemical;
		using namespace core::graph;
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

		ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( "pre_talaris_2013_standard" );
		(*sfxn)( *trpcage ); // score the pose first;
		sfxn->setup_for_packing( *trpcage, task->repacking_residues(), task->designing_residues() );
		GraphOP packer_neighbor_graph = create_packer_graph( *trpcage, *sfxn, task );

		RotamerSetsOP rotsets( new RotamerSets() );
		rotsets->set_task( task );
		rotsets->build_rotamers( *trpcage, *sfxn, packer_neighbor_graph );
		rotsets->prepare_sets_for_packing( *trpcage, *sfxn );


		FASTERInteractionGraphOP faster_ig = new FASTERInteractionGraph( 3 );
		//core::pack::pack_rotamers_setup( *trpcage, *sfxn, task, rot_sets, ig );

		rotsets->compute_energies( *trpcage, *sfxn, packer_neighbor_graph, faster_ig );

		/*for ( Size ii = 1; ii <= rotsets->nmoltenres(); ++ii ) {
			std::cout << "Rotset " << ii << " with " << rotsets->rotamer_set_for_moltenresidue(ii)->num_rotamers() << " rotamers" << std::endl;
		}*/

		//std::cout.precision( 8 );

		faster_ig->prepare_for_simulated_annealing();
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

		TS_ASSERT( netstate( 1 ) ==  1 );
		TS_ASSERT( netstate( 2 ) ==  6 );
		TS_ASSERT( netstate( 3 ) ==  4 );

		//for ( core::Size ii = 1; ii <= 3; ++ii ) {
		//	trpcage->replace_residue( ii + 10, *rotsets->rotamer_set_for_moltenresidue( ii )->rotamer( netstate(ii) ), false );
		//}
		//std::cout << "score for bmec" << (*sfxn)( *trpcage );

		TS_ASSERT_DELTA( faster_ig->get_energy_current_state_assignment(), -3.1775541, 1e-5 );

		faster_ig->relax_in_current_context();
		faster_ig->commit_relaxation();
		faster_ig->get_current_network_state( netstate );

		/*std::cout << "BMEC relaxed relaxed state:";
		for ( Size ii = 1; ii <= 3; ++ii ) {
			std::cout << " " << netstate( ii );
		}
		std::cout << std::endl;
		std::cout << faster_ig->get_energy_current_state_assignment() << std::endl;*/

		TS_ASSERT( netstate( 1 ) ==  1 );
		TS_ASSERT( netstate( 2 ) ==  1 );
		TS_ASSERT( netstate( 3 ) ==  1 );

		TS_ASSERT_DELTA( faster_ig->get_energy_current_state_assignment(), -2.8373156, 1e-5 );

		faster_ig->perturb_sBR_and_relax( 2, 6 );
		faster_ig->commit_relaxation();
		faster_ig->get_current_network_state( netstate );

		/*std::cout << "sPBR (2,6) relaxed state:";
		for ( Size ii = 1; ii <= 3; ++ii ) {
			std::cout << " " << netstate( ii );
		}
		std::cout << std::endl;
		std::cout << faster_ig->get_energy_current_state_assignment() << std::endl;*/

		TS_ASSERT( netstate( 1 ) ==  1 );
		TS_ASSERT( netstate( 2 ) ==  6 );
		TS_ASSERT( netstate( 3 ) ==  4 );

		TS_ASSERT_DELTA( faster_ig->get_energy_current_state_assignment(), -3.1775541, 1e-5 );

		core::PackerEnergy delta1 = faster_ig->perturb_sBR_and_relax( 3, 5 );
		faster_ig->commit_relaxation();
		faster_ig->get_current_network_state( netstate );

		TS_ASSERT_DELTA( faster_ig->get_energy_current_state_assignment(), -2.8532233, 1e-5 );
		TS_ASSERT_DELTA( delta1, -2.8532233 - -3.1775541, 1e-5 );

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

		TS_ASSERT( netstate( 1 ) ==  1 );
		TS_ASSERT( netstate( 2 ) ==  35 );
		TS_ASSERT( netstate( 3 ) ==  4 );

		TS_ASSERT_DELTA( faster_ig->get_energy_current_state_assignment(), -1.3659014, 1e-5 );

		core::PackerEnergy delta3, dummy;
		faster_ig->consider_substitution( 2, 6, delta3, dummy );
		TS_ASSERT_DELTA( delta3, -3.1775541 - -1.3659014, 1e-5 );

		core::PackerEnergy delta4 = faster_ig->perturb_sBR_and_relax( 3, 20 );
		faster_ig->commit_relaxation();
		faster_ig->get_current_network_state( netstate );

		/*std::cout << "sPBR (3,20) relaxed state:";
		for ( Size ii = 1; ii <= 3; ++ii ) {
			std::cout << " " << netstate( ii );
		}
		std::cout << std::endl;
		std::cout << faster_ig->get_energy_current_state_assignment() << std::endl;*/

		TS_ASSERT( netstate( 1 ) ==  1 );
		TS_ASSERT( netstate( 2 ) ==  6 );
		TS_ASSERT( netstate( 3 ) ==  20 );

		TS_ASSERT_DELTA( delta4, -1.4442203 - -1.3659014, 1e-5 );
		TS_ASSERT_DELTA( faster_ig->get_energy_current_state_assignment(), -1.4442203, 1e-5 );


	}


};
