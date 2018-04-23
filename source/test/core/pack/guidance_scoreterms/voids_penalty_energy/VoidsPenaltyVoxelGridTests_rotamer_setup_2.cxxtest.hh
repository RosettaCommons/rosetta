// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/pack/voids_penalty_energy/VoidsPenaltyVoxelGridTests.cxxtest.hh
/// @brief  Unit tests for the VoidsPenaltyVoxelGrid class.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <core/pack/guidance_scoreterms/voids_penalty_energy/VoidsPenaltyVoxelGrid.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/scoring/ScoreFunction.hh>

// Protocols Headers (for convenience, for pose manipulation).
//#include <protocols/simple_moves/MutateResidue.hh>

// Utility, etc Headers
#include <numeric/xyzVector.hh>
#include <basic/Tracer.hh>
#include <map>

static basic::Tracer TR("VoidsPenaltyVoxelGridTests_rotamer_setup_2");


class VoidsPenaltyVoxelGridTests_rotamer_setup_2 : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	void test_rotamer_setup_2() {
		core::pack::guidance_scoreterms::voids_penalty_energy::VoidsPenaltyVoxelGrid voxelgrid;
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/pack/guidance_scoreterms/voids_penalty_energy/1ubq.pdb", false, core::import_pose::PDB_file );
		voxelgrid.set_voxel_size_and_padding(0.5, 1.0);
		voxelgrid.set_up_voxel_grid_and_compute_burial(pose);

		core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( pose ));

		utility::vector1< bool > keep_aas( core::chemical::num_canonical_aas, false );
		keep_aas[ core::chemical::aa_ala ] = true;
		keep_aas[ core::chemical::aa_phe ] = true;
		keep_aas[ core::chemical::aa_met ] = true;
		keep_aas[ core::chemical::aa_ile ] = true;
		keep_aas[ core::chemical::aa_leu ] = true;
		keep_aas[ core::chemical::aa_tyr ] = true;
		keep_aas[ core::chemical::aa_val ] = true;
		keep_aas[ core::chemical::aa_trp ] = true;
		keep_aas[ core::chemical::aa_gly ] = true;

		for ( core::Size i = 1; i <= pose.total_residue(); i++ ) {
			task->nonconst_residue_task( i ).include_current();
			task->nonconst_residue_task( i ).restrict_absent_canonical_aas( keep_aas );
		}

		core::pack::rotamer_set::RotamerSetsOP rotsets( new core::pack::rotamer_set::RotamerSets() );
		rotsets->set_task( task );
		core::scoring::ScoreFunction scorefxn;
		scorefxn(pose);
		scorefxn.set_weight( core::scoring::fa_rep, 0.1 );
		utility::graph::GraphOP packer_neighbor_graph( core::pack::create_packer_graph( pose, scorefxn, task ) );
		rotsets->build_rotamers( pose, scorefxn, packer_neighbor_graph );

		utility::vector1< std::map< core::conformation::ResidueCOP, core::Real > > rotamer_volume_maps;

		voxelgrid.prune_voxels_for_fixed_residues( pose, *rotsets, nullptr, nullptr );
		voxelgrid.compute_volumes_of_buried_rotamers( pose, *rotsets, rotamer_volume_maps, nullptr, nullptr );
		TR << "Total volume is " << voxelgrid.total_buried_volume() << " A^3.  Reachable volume is " << voxelgrid.reachable_buried_volume() << " A^3." << std::endl;

		TS_ASSERT_DELTA( voxelgrid.total_buried_volume(), voxelgrid.reachable_buried_volume(), voxelgrid.total_buried_volume_*0.02 );
	}
};
