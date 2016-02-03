// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/RotamerConstraint.cxxtest.hh
/// @brief  test suite for Rotamer constraints
/// @author Rocco Moretti

// Test headers
#include <core/pack/dunbrack/RotamerConstraint.hh>

#include <cxxtest/TestSuite.h>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/io/pdb/build_pose_as_is.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <test/core/init_util.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR("core.scoring.constraints.RotamerConstraint.cxxtest");

class RotamerConstraintDun2002Tests : public CxxTest::TestSuite
{

private:
	std::string filename_; // Has on-rotamer values
	core::scoring::ScoreFunctionOP scorefxn_;

public:

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-restore_pre_talaris_2013_behavior -override_rsd_type_limit" );

		filename_ = "core/pack/dunbrack/1UBQ_repack.pdb";

		scorefxn_ = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction );
		scorefxn_->set_weight(core::scoring::fa_dun, 1.0);
		scorefxn_->set_weight(core::scoring::dunbrack_constraint, 0.000001); // Needed to turn on constraints

	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_unbound_rot() {
		using namespace core::pack::dunbrack;

		core::pose::Pose pose;
		core::io::pdb::build_pose_from_pdb_as_is(pose, filename_);
		core::Real original_score = scorefxn_->score(pose);

		core::pose::Pose unboundrot_pose(pose);
		core::pose::PoseCOPs poselist( 1, core::pose::PoseOP( new core::pose::Pose(pose) ));
		load_unboundrot(unboundrot_pose, poselist);
		core::Real unboundrot_score = scorefxn_->score(unboundrot_pose);

		// Check that repacking the original pose results in comparable (fa_dun) energy to the unboundrot values.
		core::pack::task::PackerTaskOP task(core::pack::task::TaskFactory::create_packer_task(pose));
		task->restrict_to_repacking();
		core::pack::pack_rotamers(pose, *scorefxn_, task);
		core::Real packed_score = scorefxn_->score(pose);

		TR << "Totals: " << " Original: " << original_score <<" With unboundrot: " << unboundrot_score << " Original repacked: " << packed_score << std::endl;

		TS_ASSERT_DELTA(original_score, 73.5471, 0.1);
		TS_ASSERT_DELTA(unboundrot_score, 53.9043, 0.1);
		TS_ASSERT_DELTA(unboundrot_score, packed_score, 0.1);
	}


	/// Test for a bug where the unbound rot corrections weren't being correctly applied in subsequent applications (specifically for the end residues).
	void test_rescore_and_copy() {
		using namespace core::pack::dunbrack;

		core::pose::Pose pose;
		core::io::pdb::build_pose_from_pdb_as_is(pose, filename_);

		core::pose::PoseCOPs poselist( 1, core::pose::PoseOP( new core::pose::Pose(pose) ));
		load_unboundrot(pose, poselist);

		core::Real score1 = scorefxn_->score(pose);
		core::pose::Pose pose_copy(pose);
		core::Real copy_score = scorefxn_->score(pose_copy);
		core::Real score2 = scorefxn_->score(pose);
		core::Real copy_score2 = scorefxn_->score(pose_copy);

		TR << "Score1: " << score1 << " Score2: " << score2 << " Copy1: " << copy_score << " Copy2: " << copy_score2 << std::endl;

		TS_ASSERT_DELTA(score1, copy_score, 0.01);
		TS_ASSERT_DELTA(score1, score2, 0.01);
		TS_ASSERT_DELTA(copy_score, copy_score2, 0.01);
	}

};
