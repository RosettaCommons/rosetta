// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file test/core/scoring/Ramachandran_power.cxxtest.hh
/// @brief Unit tests for the Ramachandran score term when used with the rama_power option.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <core/conformation/ppo_torsion_bin.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoringManager.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/chemical/AA.hh>
#include <basic/Tracer.hh>

using namespace std;

using core::Size;
using core::Real;
using core::pose::Pose;
using core::chemical::AA;
using core::scoring::Ramachandran;

static THREAD_LOCAL basic::Tracer TR("core.scoring.Ramachandran_power.cxxtest");

class RamachandranPowerTest : public CxxTest::TestSuite {

public:

	// Create a mock pose and fill it with artificial phi/psi values.  This pose
	// is used by several of the tests to compare against.
	void setUp() {
		core_init_with_additional_options("-rama_power 2.0");
		pose_ = core::pose::PoseOP( new core::pose::Pose );
		core::pose::make_pose_from_sequence(*pose_, "VVAV", "fa_standard", false);

		pose_->set_phi(1, 296); pose_->set_psi(1, 319);  // alpha helix
		pose_->set_phi(2, 235); pose_->set_psi(2, 138);  // beta strand
		pose_->set_phi(3,  55); pose_->set_psi(3,  42);  // left-handed helix
		pose_->set_phi(4, 125); pose_->set_psi(4, 100);  // forbidden
	}

	void tearDown() {
		pose_.reset();
	}

	// Test the score function on a handful of points, as defined in setup().
	void test_eval_rama_score_residue() {
		Ramachandran const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran();
		Real expected[] = { -0.2578, -0.9390, 1.1552, 34.6217 };

		TR << "Res\tExpected\tObserved" << std::endl;
		for ( Size i = 1; i <= pose_->total_residue(); i++ ) {
			Real observed = rama.eval_rama_score_residue(pose_->residue(i));
			TR << i << "\t" << expected[i-1] << "\t" << observed << std::endl;
			TS_ASSERT_DELTA(observed, expected[i-1], 1e-4);
		}
		TR.flush();
	}

private:
	core::pose::PoseOP pose_;

};
