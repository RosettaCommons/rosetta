// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pdb1ubq.hh>

// Package Headers
#include <protocols/simple_moves/ConcatenatePosesMover.fwd.hh>
#include <protocols/simple_moves/ConcatenatePosesMover.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/types.hh>

// Utility Headers
#include <numeric/xyzVector.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

using namespace core;

static basic::Tracer TR("protocols.simple_moves.ConcatenatePoses.cxxtest");

class ConcatenatePosesTests : public CxxTest::TestSuite {

private: // variables
	core::pose::PoseOP testpose_;

public: // test functions

	// Test Setup Functions ///////////////////////////

	/// @brief Setup Test
	void setUp() {

		using namespace core::pose;
		using namespace core::import_pose;

		// Initialize core & options system
		core_init();

		// Load in pose from pdb
		testpose_ = pdb1ubq5to13_poseop();
		std::string component_file;
	}

	/// @brief Tear Down Test
	void tearDown() {}

	// Test Methods /////////////////////////////////

	void test_concatenator() {
		using namespace protocols::simple_moves;
		ConcatenatePosesMoverOP test_mover = ConcatenatePosesMoverOP(new ConcatenatePosesMover());
		test_mover->set_component_file("./protocols/simple_moves/concatenate_inputs/component_file");
		test_mover->apply(*testpose_);
		TS_ASSERT_EQUALS(testpose_->size(),41);
		return;
	}


}; // MutateResidueTests
