// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/conformation/Residue.cxxtest.hh
/// @brief  test suite for core::conformation::Residue
/// @author Christopher Miles (cmiles@uw.edu)

// Test Headers
#include <cxxtest/TestSuite.h>

// Unit Headers
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>

// Project headers
#include <test/core/init_util.hh>

//Auto Headers
#include <utility/vector1.hh>


namespace {

class ResidueTest : public CxxTest::TestSuite {
public:
	core::pose::Pose pose_;

	void setUp() {
		core_init();
		core::import_pose::pose_from_pdb(pose_, "core/conformation/4gatA.pdb");
	}

	void test_isDNA() {
		unsigned dna_start = 68;
		unsigned dna_end = 93;

		for ( unsigned i = 1; i < dna_start; ++i ) {
			TS_ASSERT(!pose_.residue(i).is_DNA());
		}

		for ( unsigned i = dna_start; i <= dna_end; ++i ) {
			TS_ASSERT(pose_.residue(i).is_DNA());
		}
	}

	void test_isLigand() {
		TS_ASSERT(pose_.residue(67).is_ligand());
	}
};
}  // anonymous namespace
