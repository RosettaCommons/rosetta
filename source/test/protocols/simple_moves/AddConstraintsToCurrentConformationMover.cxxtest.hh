// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_moves/AddConstraintsToCurrentConformationMover.cxxtest.hh
/// @brief  test for AddConstraintsToCurrentConformationMover mover
/// @author Steven Lewis smlewi@gmail.com

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/rosettascripts.hh>
#include <test/util/pose_funcs.hh>

// Project Headers
#include <core/types.hh>

#include <protocols/simple_moves/AddConstraintsToCurrentConformationMover.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/func/Func.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>

static basic::Tracer TR("protocols.simple_moves.AddConstraintsToCurrentConformationMover.cxxtest.hh");

// --------------- Test Class --------------- //

class AddConstraintsToCurrentConformationMoverTests : public CxxTest::TestSuite {

private:
	core::pose::PoseOP test_dimer_pose_;
	core::scoring::ScoreFunctionOP scorefxn_;
public:

	void setUp() {
		core_init();

		test_dimer_pose_ = create_2res_1ten_2res_trp_cage_poseop(); //dimer structure
	}

	void tearDown() {
	}

	// //this was used for setting the tests up
	// void cst_shower() {
	//
	//  core::scoring::constraints::ConstraintSetCOP cstset( test_dimer_pose_->constraint_set());
	//
	//  //core::scoring::constraints::ConstraintCOPs
	//  auto csts(cstset->get_all_constraints());
	//  TR << csts.size() << std::endl;
	//  for (auto eachcst : csts){
	//   TR << eachcst->type() << std::endl; //show(TR);
	//   eachcst->get_func().show_definition(TR);
	//   TR << std::endl;
	//  }
	// }

	//each of these tests tests a similar thing
	void multi_test( core::Size const num_cst, std::string const & cst_type, std::string const & func_type) {

		core::scoring::constraints::ConstraintSetCOP cstset( test_dimer_pose_->constraint_set());

		TS_ASSERT(!cstset->is_empty());
		TS_ASSERT(cstset->has_constraints());
		TS_ASSERT(cstset->has_residue_pair_constraints());

		//core::scoring::constraints::ConstraintCOPs
		auto csts(cstset->get_all_constraints());
		TS_ASSERT_EQUALS(csts.size(), num_cst);  //how many constraints?
		for ( auto eachcst : csts ) {
			TS_ASSERT_EQUALS(cst_type, eachcst->type()); //what type of constraint?
			std::ostringstream func_def;
			eachcst->get_func().show_definition(func_def);
			//split definition into vector of strings; first element is func type
			//what type of func?
			TS_ASSERT_EQUALS(func_type, utility::split(func_def.str())[1]);
		}
		return;
	}


	//SML Dec 19 2016 - this is obviously not a complete test but it's better than zero

	void test_atompair_bounded() {
		protocols::simple_moves::AddConstraintsToCurrentConformationMover actccm;
		actccm.use_distance_cst() = true;
		actccm.use_bounded_func() = true;
		actccm.bound_width() = 8.0;
		actccm.max_distance() = 10000; //just something large, I have no idea how far apart the test pose residues are
		actccm.min_seq_sep() = 1;

		actccm.apply(*test_dimer_pose_);
		TR << "TEST 1" << std::endl;
		//test_dimer_pose_->constraint_set()->show(TR);
		//TR << std::endl;
		//cst_shower();
		multi_test(6, "AtomPair", "BOUNDED");

	}

	void test_atompair_harmonic() {
		protocols::simple_moves::AddConstraintsToCurrentConformationMover actccm;
		actccm.use_distance_cst() = true;
		actccm.use_bounded_func() = false;
		actccm.coord_dev() = 1.0;
		actccm.max_distance() = 10000; //just something large, I have no idea how far apart the test pose residues are
		actccm.min_seq_sep() = 1;

		actccm.apply(*test_dimer_pose_);
		TR << "TEST 2" << std::endl;
		//test_dimer_pose_->constraint_set()->show(TR);
		//TR << std::endl;
		//cst_shower();
		multi_test(6, "AtomPair", "SCALARWEIGHTEDFUNC");

	}

	void test_coord_bounded() {
		protocols::simple_moves::AddConstraintsToCurrentConformationMover actccm;
		actccm.use_distance_cst() = false;
		actccm.use_bounded_func() = true;
		actccm.bound_width() = 8.0;

		actccm.apply(*test_dimer_pose_);
		//TR << "TEST 3" << std::endl;
		//test_dimer_pose_->constraint_set()->show(TR);
		//TR << std::endl;
		//cst_shower();
		multi_test(4, "CoordinateConstraint", "BOUNDED");

	}

	void test_coord_harmonic() {
		protocols::simple_moves::AddConstraintsToCurrentConformationMover actccm;
		actccm.use_distance_cst() = false;
		actccm.use_bounded_func() = false;
		actccm.coord_dev() = 1.0;

		actccm.apply(*test_dimer_pose_);
		//TR << "TEST 4" << std::endl;
		//test_dimer_pose_->constraint_set()->show(TR);
		//TR << std::endl;
		//cst_shower();
		multi_test(4, "CoordinateConstraint", "HARMONIC");

	}

};
