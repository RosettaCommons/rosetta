// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/FASTERAnnealer.cxxtest.hh
/// @brief  test suite for the FASTER Annealer
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test framework headers
#include <cxxtest/TestSuite.h>

// Core Headers
#include <core/pack/annealer/FASTERAnnealer.hh>
#include <core/pack/interaction_graph/FASTERInteractionGraph.hh>

#include <core/chemical/AA.hh>

#include <core/graph/Graph.hh>

// AUTO-REMOVED #include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
// AUTO-REMOVED #include <core/pack/rotamer_set/RotamerSet.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>


// Test headers
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

//Auto Headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>


class FASTERAnnealerTests : public CxxTest::TestSuite {
public:

	void setUp() {
		core_init_with_additional_options( "-restore_pre_talaris_2013_behavior -override_rsd_type_limit" );
	}

	void test_instantiate_FASTER_annealer() {
		using namespace core::chemical;
		using namespace core::graph;
		using namespace core::pack;
		using namespace core::pack::annealer;
		using namespace core::pack::interaction_graph;
		using namespace core::pack::rotamer_set;
		using namespace core::pack::task;
		using namespace core::pose;
		using namespace core::scoring;
		using core::Size;

		PoseOP trpcage = create_trpcage_ideal_poseop();
		PackerTaskOP task = TaskFactory::create_packer_task( *trpcage );

		utility::vector1< bool > allowed_aas( num_canonical_aas, false );
		allowed_aas[ aa_ala ] = allowed_aas[ aa_tyr ] = allowed_aas[ aa_phe ] = allowed_aas[ aa_leu ] = allowed_aas[ aa_trp ] = true;

		for ( Size ii = 1; ii <= 20; ++ii ) {
			if ( ii == 3 || ii == 4 || ii == 6 || ii == 7 ) {
				task->nonconst_residue_task( ii ).restrict_absent_canonical_aas( allowed_aas );
				task->nonconst_residue_task( ii ).or_ex1( true );
				task->nonconst_residue_task( ii ).or_ex2( true );
				task->nonconst_residue_task( ii ).and_extrachi_cutoff( 1 );
			} else {
				task->nonconst_residue_task( ii ).prevent_repacking();
			}
		}

		ScoreFunctionOP sfxn = get_score_function();
		(*sfxn)( *trpcage ); // score the pose first;
		sfxn->setup_for_packing( *trpcage, task->repacking_residues(), task->designing_residues() );
		GraphOP packer_neighbor_graph = create_packer_graph( *trpcage, *sfxn, task );

		RotamerSetsOP rotsets( new RotamerSets() );
		rotsets->set_task( task );
		rotsets->build_rotamers( *trpcage, *sfxn, packer_neighbor_graph );
		rotsets->prepare_sets_for_packing( *trpcage, *sfxn );


		FASTERInteractionGraphOP faster_ig = new FASTERInteractionGraph( 4 );
		//FASTERInteractionGraphOP faster_ig = new FASTERInteractionGraph( 20 );
		//core::pack::pack_rotamers_setup( *trpcage, *sfxn, task, rot_sets, ig );

		rotsets->compute_energies( *trpcage, *sfxn, packer_neighbor_graph, faster_ig );

		/*for ( Size ii = 1; ii <= rotsets->nmoltenres(); ++ii ) {
			std::cout << "Rotset " << ii << " with " << rotsets->rotamer_set_for_moltenresidue(ii)->num_rotamers() << " rotamers" << std::endl;
		}*/

		//std::cout.precision( 8 );

		using namespace ObjexxFCL;

		FArray1D_int bestrotamer_at_seqpos( trpcage->total_residue() );
		core::PackerEnergy bestenergy( 0.0 );
		bool start_with_current = false;
		FArray1D_int current_rot_index( trpcage->total_residue(), 0 );
		bool calc_rot_freq = false;
		FArray1D< core::PackerEnergy > rot_freq( faster_ig->get_num_total_states(), 0.0 );

		FASTERAnnealer fa( bestrotamer_at_seqpos, bestenergy, start_with_current,
			faster_ig, rotsets, current_rot_index, calc_rot_freq, rot_freq );

		fa.run();

		//std::cout << "Best energy: " << bestenergy << std::endl;

		ObjexxFCL::FArray1D_int netstate( 4 );
		//ObjexxFCL::FArray1D_int netstate( 20 );
		faster_ig->get_current_network_state( netstate );
		/*std::cout << "FASTERAnnealer final state:";
		for ( Size ii = 1; ii <= 4; ++ii ) {
			std::cout << " " << netstate( ii );
		}
		std::cout << std::endl;*/

		/// FASTERAnnealer final state: 32 1 7 35

		TS_ASSERT( netstate( 1 ) ==  32 );
		TS_ASSERT( netstate( 2 ) ==  1 );
		TS_ASSERT( netstate( 3 ) ==  7 );
		TS_ASSERT( netstate( 4 ) ==  35 );

	}


};
