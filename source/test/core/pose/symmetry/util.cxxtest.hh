// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/pose/symmetry/util.cxxtest.hh
/// @brief  unit tests for symmetric pose utility functions
/// @author Matthew O'Meara

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <core/types.hh>
#include <protocols/simple_moves/VirtualRootMover.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>


#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "SymmetryUtilTests" );

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
		core::import_pose::pose_from_file(pose, "core/scoring/symmetry/test_in.pdb", core::import_pose::PDB_file);

		core::pose::Pose symm_pose( pose );
		std::string symm_def("core/scoring/symmetry/sym_def.dat");
		make_symmetric_pose(symm_pose, symm_def);

		ScoreFunctionOP scfxn( new ScoreFunction() );
		ScoreFunctionOP symm_scfxn( new SymmetricScoreFunction() );

		TS_ASSERT(!is_symmetric(pose));
		TS_ASSERT(!is_symmetric(*scfxn));
		TS_ASSERT(is_symmetric(symm_pose));
		TS_ASSERT(is_symmetric(*symm_scfxn));

		make_score_function_consistent_with_symmetric_state_of_pose(pose, symm_scfxn);
		TS_ASSERT(!is_symmetric(*symm_scfxn));

		make_score_function_consistent_with_symmetric_state_of_pose(symm_pose, scfxn);
		TS_ASSERT(is_symmetric(*scfxn));

	}

	void test_symmetrize_fold_tree()
	{
		using core::kinematics::FoldTree;
		using core::pose::Pose;

		Pose pose;
		core::import_pose::pose_from_file( pose, "core/scoring/symmetry/test_in.pdb", core::import_pose::PDB_file );
		Pose start_pose = pose;
		protocols::simple_moves::VirtualRootMover().apply( start_pose );

		// add symmetry data input file
		core::conformation::symmetry::SymmData symmdata1( pose.n_residue(), pose.num_jump() );
		std::string const symm_def1 = "core/conformation/symmetry/symm_def1.dat";
		symmdata1.read_symmetry_data_from_file(symm_def1);
		core::pose::symmetry::make_symmetric_pose( pose, symmdata1 );

		// try to symmeterize fold tree
		FoldTree symm_ft = start_pose.fold_tree();
		TS_ASSERT( symm_ft.check_fold_tree() );
		core::pose::symmetry::symmetrize_fold_tree( pose, symm_ft );
		TS_ASSERT( symm_ft.check_fold_tree() );
		TS_ASSERT_THROWS_NOTHING( pose.fold_tree( symm_ft ) );
	}

};
