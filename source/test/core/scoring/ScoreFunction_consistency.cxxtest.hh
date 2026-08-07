// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ScoreFunction_consistency.cxxtest.hh
/// @brief  unit tests for scorefunction consistency issues.
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Unit headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/types.hh>

#include <core/kinematics/MoveMap.hh>
#include <protocols/minimization_packing/MinMover.hh>

// Package headers
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>

#include <iostream>

static basic::Tracer TR("core.scoring.ScoreFunction.cxxtest");

// using declarations
using namespace core;
using namespace scoring;

class ScoreFunctionConsistency_beta_nov16 : public CxxTest::TestSuite {

public:
	typedef utility::keys::VariantKey< utility::options::OptionKey > VariantOptionKey;

	void setUp() {
		core_init_with_additional_options("-beta_nov16");
	}

	void tearDown() {}

	core::pose::PoseOP
	extended_pose_from_sequence(std::string const & seq) {
		core::pose::PoseOP pose = utility::pointer::make_shared< core::pose::Pose >();
		core::pose::make_pose_from_sequence(*pose,seq,"fa_standard",true);

		for ( core::Size ii(1); ii <= pose->total_residue(); ++ii ) {
			auto const & res = pose->residue(ii);
			if ( !res.is_protein() || res.is_peptoid() || res.is_carbohydrate() ) {
				continue;
			}
			pose->set_phi(ii, 180);
			pose->set_psi(ii, 180);
			pose->set_omega(ii, 180);
		}
		return pose;
	}

	void test_cloning_energy_graph() {
		core::Real const DELTA = 1e-5;
		core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();

		core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap() );
		mm->set_chi(true);
		mm->set_bb(true);

		protocols::minimization_packing::MinMover min(mm, scorefxn, "lbfgs_armijo_nonmonotone", 0.0001, true);
		min.max_iter(3000);

		core::pose::PoseOP pose_no_prescore = extended_pose_from_sequence("THISISMYTESTSTRCT");
		core::pose::PoseOP pose_prescore = extended_pose_from_sequence("THISISMYTESTSTRCT");

		core::Real prescore_value = (*scorefxn)(*pose_prescore);

		core::pose::PoseOP pose_np_cloned = pose_no_prescore->clone();
		core::pose::PoseOP pose_np_cleared = pose_no_prescore->clone();
		core::pose::PoseOP pose_p_cloned = pose_prescore->clone();
		core::pose::PoseOP pose_p_cleared = pose_prescore->clone();

		// Clearing the energies is a stand-in for deleting and then re-adding all the energy graph edges.
		pose_np_cleared->energies().clear_energies();
		pose_p_cleared->energies().clear_energies();

		TR << "Min No Prescore Cloned" << std::endl;
		min.apply(*pose_np_cloned);
		TR << "Min No Prescore Energies Cleared" << std::endl;
		min.apply(*pose_np_cleared);
		TR << "Min Prescore Cloned" << std::endl;
		min.apply(*pose_p_cloned);
		TR << "Min Prescore Energies Cleared" << std::endl;
		min.apply(*pose_p_cleared);
		TR << "Min No Prescore Original" << std::endl;
		min.apply(*pose_no_prescore);
		TR << "Min Prescore Original" << std::endl;
		min.apply(*pose_prescore);

		core::Real np_score = (*scorefxn)(*pose_no_prescore);
		core::Real p_score = (*scorefxn)(*pose_prescore);
		core::Real np_score_cloned = (*scorefxn)(*pose_np_cloned);
		core::Real p_score_cloned = (*scorefxn)(*pose_p_cloned);
		core::Real np_score_cleared = (*scorefxn)(*pose_np_cleared);
		core::Real p_score_cleared = (*scorefxn)(*pose_p_cleared);

		TS_ASSERT( p_score <= prescore_value );
		TS_ASSERT_DELTA( p_score, np_score, DELTA );
		TS_ASSERT_DELTA( p_score, np_score_cloned, DELTA );
		TS_ASSERT_DELTA( p_score, p_score_cloned , DELTA); // FAILING!! found (1406.4066 != 25.1176)
		TS_ASSERT_DELTA( p_score, np_score_cleared, DELTA );
		TS_ASSERT_DELTA( p_score, p_score_cleared , DELTA);

	}

};
