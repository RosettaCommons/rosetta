// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/SymmLinMemInteractionGraph.cxxtest.hh
/// @brief  test suite for the symmetric, minimalist on-the-fly interaction graph
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test framework headers
#include <cxxtest/TestSuite.h>

// Core Headers
#include <core/pack/interaction_graph/LinearMemoryInteractionGraph.hh>

#include <core/chemical/AA.hh>
#include <utility/graph/Graph.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>

#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/interaction_graph/PDInteractionGraph.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/make_symmetric_task.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>

// Test headers
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

//Auto Headers
#include <utility/vector1.hh>

using namespace core;

class LinMemInteractionGraphTests : public CxxTest::TestSuite {
public:

	void setUp() {
		core_init();
	}

	void test_linmem_ig() {
		using namespace core::chemical;
		using namespace core::conformation::symmetry;
		using namespace utility::graph;
		using namespace core::pack;
		using namespace core::pack::interaction_graph;
		using namespace core::pack::rotamer_set;
		using namespace core::pack::task;
		using namespace core::pose;
		using namespace core::scoring;
		using core::Size;
		using core::Real;

		Pose pose;
		core::import_pose::pose_from_file( pose, "core/scoring/symmetry/test_in.pdb" , core::import_pose::PDB_file);
		PackerTaskOP task = TaskFactory::create_packer_task( pose );

		utility::vector1< bool > allowed_aas( num_canonical_aas, false );
		allowed_aas[ aa_ala ] = allowed_aas[ aa_gly ] = allowed_aas[ aa_phe ] = allowed_aas[ aa_pro ] = allowed_aas[ aa_arg ] = true;

		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			if ( ii == 11 || ii == 12 || ii == 13 ) {
				task->nonconst_residue_task( ii ).restrict_absent_canonical_aas( allowed_aas );
			} else {
				task->nonconst_residue_task( ii ).prevent_repacking();
			}
		}

		ScoreFunctionOP sfxn( new core::scoring::ScoreFunction ); //get_score_function();
		sfxn->set_weight( fa_atr, 0.8 );
		sfxn->set_weight( fa_rep, 0.44 );
		sfxn->set_weight( fa_sol, 0.65 );
		sfxn->set_weight( lk_ball, 0.65 );
		methods::EnergyMethodOptionsOP emopts( new methods::EnergyMethodOptions( sfxn->energy_method_options() ) );
		emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );
		sfxn->set_energy_method_options( *emopts );

		//core::Real const initial_score = (*sfxn)( pose ); // score the pose first;
		(*sfxn)( pose ); // score the pose first;
		sfxn->setup_for_packing( pose, task->repacking_residues(), task->designing_residues() );
		GraphOP packer_neighbor_graph = create_packer_graph( pose, *sfxn, task );

		core::pack::rotamer_set::RotamerSetsOP rotsets( new core::pack::rotamer_set::RotamerSets() );
		rotsets->initialize_pose_for_rotsets_creation(pose); //fpd update Tsymm_
		rotsets->set_task( task );
		rotsets->build_rotamers( pose, *sfxn, packer_neighbor_graph );
		rotsets->prepare_sets_for_packing( pose, *sfxn );

		PDInteractionGraphOP regular_ig( new PDInteractionGraph( 3 ) );
		LinearMemoryInteractionGraphOP linmem_ig( new LinearMemoryInteractionGraph( 3 ) );
		linmem_ig->set_pose( pose );
		linmem_ig->set_score_function( *sfxn );

		rotsets->compute_energies( pose, *sfxn, packer_neighbor_graph, linmem_ig );
		rotsets->compute_energies( pose, *sfxn, packer_neighbor_graph, regular_ig );

		regular_ig->prepare_for_simulated_annealing();
		linmem_ig->prepare_for_simulated_annealing();


		// now let's make sure that the on-the-fly interaction graph is correctly calculating the
		// change in energies for a series of rotamer substitutions.
		{
			core::PackerEnergy deltaE, prevnode_energy;
			core::PackerEnergy regular_deltaE, regular_prevnode_energy;
			for ( Size ii = 11; ii <= 13; ++ii ) {
				pose.replace_residue( ii, *rotsets->rotamer_set_for_residue( ii )->rotamer( 1 ), false );
				linmem_ig->consider_substitution( ii-10, 1, deltaE, prevnode_energy );
				linmem_ig->commit_considered_substitution();
				regular_ig->consider_substitution( ii-10, 1, deltaE, prevnode_energy );
				regular_ig->commit_considered_substitution();
			}
			core::Real state1_energy = (*sfxn)( pose );
			for ( Size ii = 11; ii <= 13; ++ii ) {
				Size iimoltenresid = ii-10;
				//std::cout << "testing rotamers from residue " << ii << std::endl;
				for ( Size jj = 1; jj <= rotsets->rotamer_set_for_residue( ii )->num_rotamers(); ++jj ) {
					//std::cout << "  rotamer #" << jj << std::endl;
					pose.replace_residue( ii, *rotsets->rotamer_set_for_residue( ii )->rotamer( jj ), false );
					core::Real jjenergy = (*sfxn)( pose );

					linmem_ig->consider_substitution( iimoltenresid, jj, deltaE, prevnode_energy );
					regular_ig->consider_substitution( iimoltenresid, jj, regular_deltaE, regular_prevnode_energy );


					TS_ASSERT( std::abs( ( deltaE - (jjenergy - state1_energy)) / ( jjenergy - state1_energy + 1e-6 ) ) < 2e-4 || ( std::abs( jjenergy - state1_energy ) < 1.0 && std::abs( deltaE - ( jjenergy - state1_energy )) < 2e-4 ));
					TS_ASSERT( std::abs( ( deltaE - regular_deltaE) / ( regular_deltaE + 1e-6 ) ) < 2e-4 || ( std::abs( regular_deltaE ) < 1.0 && std::abs( deltaE - regular_deltaE ) < 2e-4 ));
				}
				// restore pose to state1 assigment
				pose.replace_residue( ii, *rotsets->rotamer_set_for_residue( ii )->rotamer( 1 ), false );
			}
		}

	}


};
