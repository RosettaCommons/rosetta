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

static basic::Tracer TR("VoidsPenaltyVoxelGridTests");


class VoidsPenaltyVoxelGridTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	/// @brief Test the VoidsPenaltyVoxelGrid::compute_bounding_box() function.
	void test_compute_bounding_box() {
		core::pack::guidance_scoreterms::voids_penalty_energy::VoidsPenaltyVoxelGrid voxelgrid;
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/pack/guidance_scoreterms/voids_penalty_energy/1ubq.pdb", false, core::import_pose::PDB_file );
		numeric::xyzVector< core::Real > lower_left_coords(0.0, 0.0, 0.0);
		numeric::xyzVector< core::Real > upper_right_coords(0.0, 0.0, 0.0);
		voxelgrid.compute_bounding_box( pose, 1.0, lower_left_coords, upper_right_coords );

		for ( core::Size ir(1), irmax(pose.total_residue()); ir<=irmax; ++ir ) {
			for ( core::Size ia(1), iamax(pose.residue_type(ir).natoms()); ia<=iamax; ++ia ) {
				numeric::xyzVector< core::Real > curpos( pose.xyz( core::id::AtomID(ia, ir) ) );
				TS_ASSERT_LESS_THAN( lower_left_coords.x(), curpos.x() );
				TS_ASSERT_LESS_THAN( lower_left_coords.y(), curpos.y() );
				TS_ASSERT_LESS_THAN( lower_left_coords.z(), curpos.z() );
				TS_ASSERT_LESS_THAN( curpos.x(), upper_right_coords.x() );
				TS_ASSERT_LESS_THAN( curpos.y(), upper_right_coords.y() );
				TS_ASSERT_LESS_THAN( curpos.z(), upper_right_coords.z() );
			}
		}
	}

	void test_voxel_grid() {
		core::pack::guidance_scoreterms::voids_penalty_energy::VoidsPenaltyVoxelGrid voxelgrid;
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/pack/guidance_scoreterms/voids_penalty_energy/1ubq.pdb", false, core::import_pose::PDB_file );
		voxelgrid.set_voxel_size_and_padding(1.5, 1.5);
		voxelgrid.set_up_voxel_grid_and_compute_burial(pose);

		TS_ASSERT_DELTA( voxelgrid.lower_left_coords_.x(), 13.94, 1e-2 );
		TS_ASSERT_DELTA( voxelgrid.lower_left_coords_.y(), 11.3861, 1e-2 );
		TS_ASSERT_DELTA( voxelgrid.lower_left_coords_.z(), -1.866, 1e-2 );
		TS_ASSERT_DELTA( voxelgrid.upper_right_coords_.x(), 48.44, 1e-2 );
		TS_ASSERT_DELTA( voxelgrid.upper_right_coords_.y(), 47.3861, 1e-2 );
		TS_ASSERT_DELTA( voxelgrid.upper_right_coords_.z(), 38.6340, 1e-2 );

		TS_ASSERT_EQUALS( voxelgrid.voxel_data_dimensions_[1], 24 );
		TS_ASSERT_EQUALS( voxelgrid.voxel_data_dimensions_[2], 25 );
		TS_ASSERT_EQUALS( voxelgrid.voxel_data_dimensions_[3], 28 );

		//TS_ASSERT_EQUALS( voxelgrid.visualize_voxel_grid(false)->total_residue(), 14904 );
		TS_ASSERT_EQUALS( voxelgrid.voxel_data_dimensions_[1]*voxelgrid.voxel_data_dimensions_[2]*voxelgrid.voxel_data_dimensions_[3], 16800 );
		//TS_ASSERT_EQUALS( voxelgrid.visualize_voxel_grid(true)->total_residue(), 589);
		TS_ASSERT_DELTA( voxelgrid.voxel_size_, 1.5, 1e-4);
		TS_ASSERT_DELTA( voxelgrid.half_voxel_size_, 0.75, 1e-4);
		TS_ASSERT_DELTA( voxelgrid.voxel_volume_, 3.375, 1e-4);
		TS_ASSERT_DELTA( voxelgrid.total_buried_volume(), 555.0*pow(1.5, 3), 1e-4 );

		numeric::xyzVector< core::Real > const testpoint1( 13.4, 11.2, -1.9 ); //Should be in [0,0,0]
		numeric::xyzVector< core::Real > const testpoint2( 15.64, 11.2, -1.9 ); //Should be in [1,0,0]
		numeric::xyzVector< core::Real > const testpoint3( 48.3, 47.4, 37.9 ); //Should be in [max, max, max]=[23,24,27].
		numeric::xyzVector< core::Real > const testpoint4( 3.0, 11.2, -1.9 ); //Should be outside box.
		numeric::xyzVector< core::Real > const testpoint5( 15.64, 1.2, -1.9 ); //Should be outside box.
		numeric::xyzVector< core::Real > const testpoint6( 13.0, 12.7, -10.9 ); //Should be outside box.
		numeric::xyzVector< core::Real > const testpoint7( 57.8, 47.4, 37.1 ); //Should be outside box.
		numeric::xyzVector< core::Real > const testpoint8( 46.9, 56.8, 38.5 ); //Should be outside box.
		numeric::xyzVector< core::Real > const testpoint9( 48.3, 45.9, 47.9 ); //Should be outside box.

		utility::fixedsizearray1< core::Size, 3 > indices1(1111);
		utility::fixedsizearray1< core::Size, 3 > indices2(2222);
		utility::fixedsizearray1< core::Size, 3 > indices3(3333);
		utility::fixedsizearray1< core::Size, 3 > indices4(4444);
		utility::fixedsizearray1< core::Size, 3 > indices5(5555);
		utility::fixedsizearray1< core::Size, 3 > indices6(6666);
		utility::fixedsizearray1< core::Size, 3 > indices7(7777);
		utility::fixedsizearray1< core::Size, 3 > indices8(8888);
		utility::fixedsizearray1< core::Size, 3 > indices9(9999);

		TS_ASSERT( voxelgrid.get_indices_of_voxel_from_coordinates( testpoint1, indices1 ) ); //True means inside box.
		TS_ASSERT( voxelgrid.get_indices_of_voxel_from_coordinates( testpoint2, indices2 ) ); //True means inside box.
		TS_ASSERT( voxelgrid.get_indices_of_voxel_from_coordinates( testpoint3, indices3 ) ); //True means inside box.
		TS_ASSERT( !voxelgrid.get_indices_of_voxel_from_coordinates( testpoint4, indices4 ) ); //False means outside box.
		TS_ASSERT( !voxelgrid.get_indices_of_voxel_from_coordinates( testpoint5, indices5 ) ); //False means outside box.
		TS_ASSERT( !voxelgrid.get_indices_of_voxel_from_coordinates( testpoint6, indices6 ) ); //False means outside box.
		TS_ASSERT( !voxelgrid.get_indices_of_voxel_from_coordinates( testpoint7, indices7 ) ); //False means outside box.
		TS_ASSERT( !voxelgrid.get_indices_of_voxel_from_coordinates( testpoint8, indices8 ) ); //False means outside box.
		TS_ASSERT( !voxelgrid.get_indices_of_voxel_from_coordinates( testpoint9, indices9 ) ); //False means outside box.

		TS_ASSERT_EQUALS( indices1[1], 0 ); TS_ASSERT_EQUALS( indices1[2], 0 ); TS_ASSERT_EQUALS( indices1[3], 0 );
		TS_ASSERT_EQUALS( indices2[1], 1 ); TS_ASSERT_EQUALS( indices2[2], 0 ); TS_ASSERT_EQUALS( indices2[3], 0 );
		TS_ASSERT_EQUALS( indices3[1], 23 ); TS_ASSERT_EQUALS( indices3[2], 24 ); TS_ASSERT_EQUALS( indices3[3], 27 );
		TS_ASSERT_EQUALS( indices4[1], 0 ); TS_ASSERT_EQUALS( indices4[2], 0 ); TS_ASSERT_EQUALS( indices4[3], 0 );
		TS_ASSERT_EQUALS( indices5[1], 1 ); TS_ASSERT_EQUALS( indices5[2], 0 ); TS_ASSERT_EQUALS( indices5[3], 0 );
		TS_ASSERT_EQUALS( indices6[1], 0 ); TS_ASSERT_EQUALS( indices6[2], 1 ); TS_ASSERT_EQUALS( indices6[3], 0 );
		TS_ASSERT_EQUALS( indices7[1], 23 ); TS_ASSERT_EQUALS( indices7[2], 24 ); TS_ASSERT_EQUALS( indices7[3], 26 );
		TS_ASSERT_EQUALS( indices8[1], 22 ); TS_ASSERT_EQUALS( indices8[2], 24 ); TS_ASSERT_EQUALS( indices8[3], 27 );
		TS_ASSERT_EQUALS( indices9[1], 23 ); TS_ASSERT_EQUALS( indices9[2], 23 ); TS_ASSERT_EQUALS( indices9[3], 27 );
	}

};
