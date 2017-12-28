// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/pack/voids_penalty_energy/VoidsPenaltyEnergyTests_symm.cxxtest.hh
/// @brief  Unit tests for the VoidsPenaltyEnergy, a score term penalizing buried holes during packing.  This
/// suite tests the symmetric case.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <core/pack/voids_penalty_energy/VoidsPenaltyEnergy.hh>

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
#include <core/pack/make_symmetric_task.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>

// Protocols Headers -- for convenience in setting up test case.
#include <protocols/helical_bundle/MakeBundle.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("VoidsPenaltyEnergyTests");


class VoidsPenaltyEnergyTests_symm : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options( "-voids_penalty_energy_voxel_size 0.5 -voids_penalty_energy_disabled_except_during_packing false" );

	}

	void tearDown(){

	}


	/// @brief Confirm that packing with the voids_penalty score term results in low voids_penalty energy values.
	void test_packing_success(){
		core::pose::Pose pose;
		protocols::helical_bundle::MakeBundle makebundle;
		makebundle.set_default_helix_length(25);
		makebundle.set_default_r0( 7.25 );
		makebundle.set_default_omega0( -0.03490659 );
		makebundle.set_default_delta_omega1_all(0.78539816);
		makebundle.set_default_delta_omega0(0.7853981635);
		makebundle.add_helix();
		makebundle.apply(pose);

		protocols::simple_moves::symmetry::SetupForSymmetryMover set_up_d2( "core/pack/voids_penalty_energy/d2.symm" );
		set_up_d2.apply(pose);

		core::select::residue_selector::LayerSelectorOP layer_selector( new core::select::residue_selector::LayerSelector );
		layer_selector->set_use_sc_neighbors(true);
		layer_selector->set_layers( true, true, false );
		utility::vector1< bool > selection( layer_selector->apply( pose ) );

		core::scoring::symmetry::SymmetricScoreFunction sfxn;
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
		keep_aas[ core::chemical::aa_met ] = true;
		keep_aas[ core::chemical::aa_ile ] = true;
		keep_aas[ core::chemical::aa_trp ] = true;

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

		task->request_symmetrize_by_intersection();
		core::pack::task::PackerTaskOP symmtask( core::pack::make_new_symmetric_PackerTask_by_requested_method( pose, task ) );

		core::pack::rotamer_set::RotamerSetsOP rotsets( new core::pack::rotamer_set::RotamerSets() );
		rotsets->set_task( symmtask );
		utility::graph::GraphOP packer_neighbor_graph( core::pack::create_packer_graph( pose, sfxn, symmtask ) );
		rotsets->build_rotamers( pose, sfxn, packer_neighbor_graph );

		core::pose::Pose pose_copy( pose );
		sfxn(pose_copy);
		core::pack::pack_rotamers( pose_copy, sfxn, symmtask );
		post_pack_score = sfxn(pose_copy);
		TR << "Final score: " << post_pack_score << std::endl;
		TS_ASSERT_LESS_THAN( post_pack_score, 160 );

		//DELETE THE FOLLOWING: for debugging only!
		pose.dump_pdb("voids_symm_prepack.pdb");
		pose_copy.dump_pdb("voids_symm_postpack.pdb");

	}



};
