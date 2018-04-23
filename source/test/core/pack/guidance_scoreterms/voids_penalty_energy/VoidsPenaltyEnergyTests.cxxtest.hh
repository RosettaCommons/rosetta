// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/pack/voids_penalty_energy/VoidsPenaltyEnergyTests.cxxtest.hh
/// @brief  Unit tests for the VoidsPenaltyEnergy, a score term penalizing buried holes during packing.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <core/pack/guidance_scoreterms/voids_penalty_energy/VoidsPenaltyEnergy.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/LayerSelector.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>


// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("VoidsPenaltyEnergyTests");


class VoidsPenaltyEnergyTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options( "-voids_penalty_energy_voxel_size 0.5" );

	}

	void tearDown(){

	}


	/// @brief Confirm that packing with the voids_penalty score term results in low voids_penalty energy values.
	void test_packing_success(){
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/pack/guidance_scoreterms/voids_penalty_energy/1ubq.pdb", false, core::import_pose::PDB_file );

		core::select::residue_selector::LayerSelectorOP layer_selector( new core::select::residue_selector::LayerSelector );
		layer_selector->set_use_sc_neighbors(true);
		layer_selector->set_layers( true, true, false );
		utility::vector1< bool > selection( layer_selector->apply( pose ) );

		core::scoring::methods::EnergyMethodOptions opts;
		opts.initialize_from_options();
		opts.voids_penalty_energy_disabled_except_during_packing(false);
		core::scoring::ScoreFunction sfxn;
		sfxn.set_energy_method_options( opts );
		//sfxn.set_weight( core::scoring::fa_rep, 0.05 );
		sfxn.set_weight( core::scoring::voids_penalty, 1.0 );
		core::Real const initial_score( sfxn(pose) );
		TR << "Intial score: " << initial_score << std::endl;
		TS_ASSERT_LESS_THAN( 2000.0, initial_score );

		utility::vector1< bool > keep_aas( core::chemical::num_canonical_aas, false );
		keep_aas[ core::chemical::aa_ala ] = true;
		keep_aas[ core::chemical::aa_gly ] = true;
		keep_aas[ core::chemical::aa_leu ] = true;
		keep_aas[ core::chemical::aa_val ] = true;
		keep_aas[ core::chemical::aa_phe ] = true;

		core::Real post_pack_score;

		//for(core::Size j(1); j<=5; ++j) {
		core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( pose ));

		for ( core::Size i = 1; i <= pose.size(); i++ ) {
			if ( selection[i] ) {
				task->nonconst_residue_task( i ).restrict_absent_canonical_aas( keep_aas );
			} else {
				task->nonconst_residue_task( i ).prevent_repacking();
			}
		}

		core::pack::rotamer_set::RotamerSetsOP rotsets( new core::pack::rotamer_set::RotamerSets() );
		rotsets->set_task( task );
		utility::graph::GraphOP packer_neighbor_graph( core::pack::create_packer_graph( pose, sfxn, task ) );
		rotsets->build_rotamers( pose, sfxn, packer_neighbor_graph );

		core::pose::Pose pose_copy( pose );
		sfxn(pose_copy);
		core::pack::pack_rotamers( pose_copy, sfxn, task );
		post_pack_score = sfxn(pose_copy);
		TR << "Final score: " << post_pack_score << std::endl;
		TS_ASSERT_LESS_THAN( post_pack_score, 160 );
		//}

	}



};
