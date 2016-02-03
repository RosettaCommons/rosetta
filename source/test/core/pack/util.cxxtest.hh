// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/util.cxxtest.hh
/// @brief  Test suite for simple util functions
/// @author JKLeman (julia.koehler1982@gmail.com)

// Package headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/util.hh>

// Test headers
#include <cxxtest/TestSuite.h>

#include <platform/types.hh>
#include <core/types.hh>

#include <test/core/init_util.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <utility/vector1.hh>

using namespace core;
using namespace core::pack;
using namespace core::pose;

class RotamerUtils : public CxxTest::TestSuite {

public:

	void setUp() {

		core_init();
	}

	void tearDown() {}


	///////////////////////////////////////////////////////////////////////////////
	/// @brief test simple rotamer recovery
	void test_rotamer_recovery() {

		TS_TRACE("Testing rotamer recovery");

		// get poses
		Pose pose, native;
		core::import_pose::pose_from_file( pose, "core/pack/4MRS__opm_0001.pdb" , core::import_pose::PDB_file);
		core::import_pose::pose_from_file( native, "core/pack/4MRS__opm.pdb" , core::import_pose::PDB_file);

		// compute rotamer recovery
		core::Real difference = 10.0;
		core::Real rot_rec( core::pack::rotamer_recovery( pose, native, difference ) );

		// test
		TS_ASSERT_DELTA( rot_rec, 0.757092, 0.00001 );

	}

	/// @brief test residue rotamer recovery
	void test_residue_rotamer_recovery() {

		TS_TRACE("Testing residue rotamer recovery");

		// get poses
		Pose pose, native;
		core::import_pose::pose_from_file( pose, "core/pack/4MRS__opm_0001.pdb" , core::import_pose::PDB_file);
		core::import_pose::pose_from_file( native, "core/pack/4MRS__opm.pdb" , core::import_pose::PDB_file);

		// compute residue rotamer recovery
		core::Real difference = 10.0;
		core::Real res_rot_rec( core::pack::residue_rotamer_recovery( pose, native, difference ) );

		// test
		TS_ASSERT_DELTA( res_rot_rec, 0.661937, 0.00001 );
	}


};
