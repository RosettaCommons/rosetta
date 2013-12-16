// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// Copyright in the Rosetta software belongs to the developers and their institutions.
// For more information, see www.rosettacommons.org.
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/pose/symmetry/util.cxxtest.hh
/// @brief  unit tests for symmetric pose utility functions
/// @author Matthew O'Meara

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>


class SymmetryUtilTests : public CxxTest::TestSuite
{


public: //setup


	SymmetryUtilTests() {};


	// Shared initialization.
	void setUp() {
		core_init();
	}


	// Shared finalization.
	void tearDown() {
	}


public: // tests

	/// @brief test Pose observer interface
	void test_make_score_function_consistent() {
		using namespace core::pose;
		using namespace core::pose::symmetry;
		using namespace core::scoring;
		using namespace core::scoring::symmetry;

		Pose pose;
		core::import_pose::pose_from_pdb(pose, "core/scoring/symmetry/test_in.pdb");

		core::pose::Pose symm_pose( pose );
		std::string symm_def("core/scoring/symmetry/sym_def.dat");
		make_symmetric_pose(symm_pose, symm_def);

		ScoreFunctionOP scfxn = new ScoreFunction();
		ScoreFunctionOP symm_scfxn = new SymmetricScoreFunction();

		TS_ASSERT(!is_symmetric(pose));
		TS_ASSERT(!is_symmetric(*scfxn));
		TS_ASSERT(is_symmetric(symm_pose));
		TS_ASSERT(is_symmetric(*symm_scfxn));

		make_score_function_consistent_with_symmetric_state_of_pose(pose, symm_scfxn);
		TS_ASSERT(!is_symmetric(*symm_scfxn));

		make_score_function_consistent_with_symmetric_state_of_pose(symm_pose, scfxn);
		TS_ASSERT(is_symmetric(*scfxn));

	}
};
