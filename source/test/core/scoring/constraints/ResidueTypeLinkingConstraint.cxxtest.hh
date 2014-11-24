// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constraints/ResidueTypeLinkingConstraint.cxxtest.hh
/// @brief  test suite for residue type linking constraints
/// @author Alex Sevy

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>
#include <test/util/deriv_funcs.hh>

#include <core/scoring/constraints/ResidueTypeLinkingConstraint.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>


class ResidueTypeLinkingConstraintTests : public CxxTest::TestSuite
{

public:
	ResidueTypeLinkingConstraintTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_constraints()
	{
		using namespace core;
		pose::Pose start_pose( create_test_in_pdb_pose() );
		pose::Pose pose1 (start_pose);
		scoring::ScoreFunctionOP scorefxn = new scoring::ScoreFunction;
		scorefxn->set_weight( scoring::res_type_linking_constraint, 1.0);
		pose1.add_constraint(new scoring::constraints::ResidueTypeLinkingConstraint(
				pose1, 1, 10, 1.0
				));
		( *scorefxn )( pose1 );
		TS_ASSERT_EQUALS( pose1.energies().total_energies()[ scoring::res_type_linking_constraint ], 0 );
		pose1.add_constraint(new scoring::constraints::ResidueTypeLinkingConstraint(
				pose1, 1, 2, 1.0
				));
		( *scorefxn )( pose1 );
		TS_ASSERT_EQUALS( pose1.energies().total_energies()[ scoring::res_type_linking_constraint ], 1 );
		pose1.add_constraint(new scoring::constraints::ResidueTypeLinkingConstraint(
				pose1, 1, 3, 2.0
				));
		( *scorefxn )( pose1 );
		TS_ASSERT_EQUALS( pose1.energies().total_energies()[ scoring::res_type_linking_constraint ], 3 );
	}


};
