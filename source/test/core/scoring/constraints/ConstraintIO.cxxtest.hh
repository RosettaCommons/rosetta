// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constaints/ConstraintIO.cxxtest.hh
/// @brief  test suite for ConstraintIO
/// @author Daniel Farrell (danpf@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>

// #include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>
#include <test/util/pdb1rpb.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>

// basic headers
#include <basic/Tracer.hh>

// numeric headers
#include <numeric/conversions.hh>

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif

using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.constraints.ConstraintIO.cxxtest.hh");

using namespace core;

class ConstraintIOTests : public CxxTest::TestSuite
{

public:
	ConstraintIOTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_constraint_read_in() {
		core::pose::Pose pose(pdb1rpb_pose());
		{
			std::string const simple_cst_txt("\n"
				"AtomPair CB 1 CB 2 SCALARWEIGHTEDFUNC 5  SUMFUNC 2 SIGMOID 35.0 0.6 CONSTANTFUNC -0.5\n"
				"AtomPair CB 4 CA 5 SCALARWEIGHTEDFUNC 5  SUMFUNC 2 SIGMOID 35.0 0.6 CONSTANTFUNC -0.5\n"
				"AtomPair CB 7 CB 10 SCALARWEIGHTEDFUNC 5  SUMFUNC 2 SIGMOID 35.0 0.6 CONSTANTFUNC -0.5\n\n"
				"AtomPair CB 8 CB 11 SCALARWEIGHTEDFUNC 5  SUMFUNC 2 SIGMOID 35.0 0.6 CONSTANTFUNC -0.5\n");
			std::stringstream ss;
			ss << simple_cst_txt;
			core::scoring::constraints::ConstraintSetOP constraint_set(
				core::scoring::constraints::ConstraintIO::get_instance()->read_constraints(
				ss,
				utility::pointer::make_shared< core::scoring::constraints::ConstraintSet >(),
				pose
				)
			);
			TS_ASSERT_EQUALS(constraint_set->get_all_constraints().size(), 4);
		}
		{
			std::string const simple_cst_txt(""
				"AtomPair CB 1 CB 2 SCALARWEIGHTEDFUNC 5  SUMFUNC 2 SIGMOID 35.0 0.6 CONSTANTFUNC -0.5\n"
				"AtomPair CB 4 CA 5 SCALARWEIGHTEDFUNC 5  SUMFUNC 2 SIGMOID 35.0 0.6 CONSTANTFUNC -0.5\n"
				"AtomPair CB 7 CB 10 SCALARWEIGHTEDFUNC 5  SUMFUNC 2 SIGMOID 35.0 0.6 CONSTANTFUNC -0.5\n\n"
				"AtomPair CB 8 CB 11 SCALARWEIGHTEDFUNC 5  SUMFUNC 2 SIGMOID 35.0 0.6 CONSTANTFUNC -0.5\n");
			std::stringstream ss;
			ss << simple_cst_txt;
			core::scoring::constraints::ConstraintSetOP constraint_set(
				core::scoring::constraints::ConstraintIO::get_instance()->read_constraints(
				ss,
				utility::pointer::make_shared< core::scoring::constraints::ConstraintSet >(),
				pose
				)
			);
			TS_ASSERT_EQUALS(constraint_set->get_all_constraints().size(), 4);
		}
	}


};
