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
#include <core/pack/voids_penalty_energy/VoidsPenaltyVoxelGrid.hh>

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

static basic::Tracer TR("VoidsPenaltyVoxelGridTests_rotamer_setup");

class VoidsPenaltyVoxelGridTests_rotamer_setup : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	void test_rotamer_setup() {
		core::pack::voids_penalty_energy::VoidsPenaltyVoxelGrid voxelgrid;
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/pack/voids_penalty_energy/1ubq.pdb", false, core::import_pose::PDB_file );
		voxelgrid.set_voxel_size_and_padding(0.5, 1.0);
		voxelgrid.set_up_voxel_grid_and_compute_burial(pose);

		core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( pose ));

		utility::vector1< bool > keep_aas( core::chemical::num_canonical_aas, false );
		keep_aas[core::chemical::aa_gly] = true;
		keep_aas[core::chemical::aa_ile] = true;

		for ( core::Size i = 1; i <= pose.total_residue(); i++ ) {
			task->nonconst_residue_task( i ).include_current();
			//task->nonconst_residue_task( i ).prevent_repacking();
			task->nonconst_residue_task( i ).restrict_absent_canonical_aas( keep_aas );
		}

		core::pack::rotamer_set::RotamerSetsOP rotsets( new core::pack::rotamer_set::RotamerSets() );
		rotsets->set_task( task );
		core::scoring::ScoreFunction scorefxn;
		scorefxn(pose);
		scorefxn.set_weight( core::scoring::fa_rep, 0.1 );
		utility::graph::GraphOP packer_neighbor_graph( core::pack::create_packer_graph( pose, scorefxn, task ) );
		rotsets->build_rotamers( pose, scorefxn, packer_neighbor_graph );

		TR << "The rotsets object contains rotamers for " << rotsets->nrotamers() << " rotamers at " << rotsets->total_residue() << " positions." << std::endl;
		TS_ASSERT_EQUALS(rotsets->total_residue(), pose.total_residue());
		TR << "Position 36 has " << rotsets->rotamer_set_for_residue(36)->num_rotamers() << " rotamers." << std::endl;

		utility::vector1< std::map< core::conformation::ResidueCOP, core::Real > > rotamer_volume_maps;
		voxelgrid.prune_voxels_for_fixed_residues(pose, *rotsets, nullptr, nullptr);
		voxelgrid.compute_volumes_of_buried_rotamers( *rotsets, rotamer_volume_maps, nullptr, nullptr );

		TS_ASSERT_EQUALS( rotamer_volume_maps.size(), pose.total_residue() );
		TS_ASSERT_EQUALS( rotamer_volume_maps[36].size(), rotsets->rotamer_set_for_residue(36)->num_rotamers() );

		bool ile_found(false), gly_found(false);

		TR << "Trying directly from map iterator:" << std::endl;
		for ( std::map< core::conformation::ResidueCOP, core::Real >::const_iterator pos36iterator( rotamer_volume_maps[36].begin() ); pos36iterator != rotamer_volume_maps[36].end(); ++pos36iterator ) {
			TR << "Volume of " << pos36iterator->first->name3() << pos36iterator->first->seqpos() << ": " << pos36iterator->second << std::endl;
			if ( pos36iterator->first->type().aa() == core::chemical::aa_gly ) {
				gly_found=true;
				TS_ASSERT_EQUALS( pos36iterator->first->name3(), "GLY" );
				TS_ASSERT_LESS_THAN( pos36iterator->second, 10.0 );
				TS_ASSERT_LESS_THAN( 0.0, pos36iterator->second );
			} else {
				ile_found=true;
				TS_ASSERT_EQUALS( pos36iterator->first->name3(), "ILE" );
				TS_ASSERT_LESS_THAN( pos36iterator->second, 65.0 );
				TS_ASSERT_LESS_THAN( 30.0, pos36iterator->second );
			}
			TS_ASSERT_EQUALS( pos36iterator->first->seqpos(), 36 );
		}
		TS_ASSERT(ile_found);
		TS_ASSERT(gly_found);
		ile_found=false;
		gly_found=false;

		TR << "Trying using RotamerCOP for map lookup:" << std::endl;
		core::pack::rotamer_set::RotamerSetCOP rotset( rotsets->rotamer_set_for_residue(36) );
		std::map< core::conformation::ResidueCOP, core::Real > const &curmap( rotamer_volume_maps[36] );

		for ( core::Size i(1), imax(rotset->num_rotamers()); i<=imax; ++i ) {
			core::conformation::ResidueCOP curres( rotset->rotamer(i) );
			TS_ASSERT_EQUALS( curmap.count( curres ), 1 );
			TR << "Volume of " << curres->name3() << curres->seqpos() << ": " << curmap.at(curres) << std::endl;
			if ( curres->aa() == core::chemical::aa_gly ) {
				gly_found=true;
				TS_ASSERT_LESS_THAN( curmap.at(curres), 10.0 );
				TS_ASSERT_LESS_THAN( 0.0, curmap.at(curres) );
			} else {
				ile_found=true;
				TS_ASSERT_EQUALS( curres->name3(), "ILE" );
				TS_ASSERT_LESS_THAN( curmap.at(curres), 65.0 );
				TS_ASSERT_LESS_THAN( 30.0, curmap.at(curres) );
			}
			TS_ASSERT_EQUALS( curres->seqpos(), 36 );
		}
		TS_ASSERT(ile_found);
		TS_ASSERT(gly_found);

		TR << "Total volume is " << voxelgrid.total_buried_volume() << " A^3.  Reachable volume is " << voxelgrid.reachable_buried_volume() << " A^3." << std::endl;
	}

};
