// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/features/CoordSerializationTests.cxxtest.hh
/// @brief test suite for serializing and deserializing xyz coordinates.  tests functions needed for compact residue schema
/// @author Sam DeLuca <samuel.l.deluca@vanderbilt.edu>

#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <protocols/features/util.hh>
#include <utility/vector1.hh>

class CoordSerializationTests : public CxxTest::TestSuite {

public:
	void setUp() {
		protocols_init();
	}

	void test_serialization() {
		core::pose::Pose pose;
		core::import_pose::pose_from_file(pose, "core/conformation/4gatA.pdb", core::import_pose::PDB_file);
		core::conformation::Residue original(pose.residue(1));
		std::string serialized_data(protocols::features::serialize_residue_xyz_coords(original));
		utility::vector1< numeric::xyzVector<core::Real> > coord_data(protocols::features::deserialize_xyz_coords(serialized_data,original.natoms()));
		TS_ASSERT_EQUALS(original.natoms(),coord_data.size());

		for ( core::Size atom_index = 1; atom_index <= original.natoms(); ++atom_index ) {
			TS_ASSERT_EQUALS(original.xyz(atom_index),coord_data[atom_index]);
		}
	}
};
