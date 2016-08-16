// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_kinematic_closure_FooBar_CXXTEST_HH
#define INCLUDED_protocols_kinematic_closure_FooBar_CXXTEST_HH

// Headers {{{1
#include <cxxtest/TestSuite.h>
#include <test/protocols/kinematic_closure/utilities.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>

// Protocol headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/utilities.hh>
#include <protocols/kinematic_closure/ClosureProblem.hh>
#include <protocols/kinematic_closure/ClosureSolution.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.hh>
#include <protocols/kinematic_closure/pivot_pickers/LoopPivots.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

// Utility headers
#include <boost/foreach.hpp>
#include <numeric/conversions.hh>
#include <basic/Tracer.hh>

#define foreach BOOST_FOREACH

// Global names {{{1
using namespace core;
using core::PointPosition;
using core::pose::Pose;
using core::id::AtomID;
using core::import_pose::pose_from_file;
using core::kinematics::FoldTree;
using protocols::kinematic_closure::ClosureProblemOP;
using protocols::kinematic_closure::ClosureProblemCOP;
using protocols::kinematic_closure::ClosureSolutionCOP;
using protocols::kinematic_closure::IndexList;
using protocols::kinematic_closure::SolutionList;
using protocols::kinematic_closure::perturbers::Perturber;
using protocols::kinematic_closure::pivot_pickers::PivotPickerOP;
using protocols::kinematic_closure::pivot_pickers::LoopPivots;
using protocols::loops::Loop;
using protocols::loops::Loops;
using numeric::conversions::DEGREES;

static THREAD_LOCAL basic::Tracer TR("protocols.kinematic_closure.TestHelpers.cxxtest");

// Forward declarations {{{1
class ClosureTest;
typedef utility::pointer::shared_ptr<ClosureTest> ClosureTestOP;
typedef utility::pointer::shared_ptr<ClosureTest const> ClosureTestCOP;
// }}}1

class ClosureTest : public Perturber { // {{{1

public:
	ClosureTest(string name, string pdb_path, Loop loop, Size num_solutions);
	~ClosureTest() {};
	string get_name() const { return name; }

public:
	virtual void perturb_subset(
			Pose const & pose,
			IndexList const & residues,
			ClosureProblemOP problem);

	virtual void perturb_for_test(ClosureProblemOP) {}

public:
	void test_closure(SolutionList) const;
	void test_closure(Pose const & pose) const;
	void test_restore(ClosureProblemCOP, SolutionList) const;
	void test_jacobian(SolutionList solutions) const;
	void test_pickers(ClosureProblemCOP, SolutionList) const;

public:
	AtomID id_from_index(Size index) const;
	Size index_from_residue(Size residue) const;

public:
	string name;
	Pose pose, original_pose;
	Loop loop;
	Size num_solutions;
	vector<Real> jacobians;
	PivotPickerOP pivot_picker;
	Real precision;
};

ClosureTest::ClosureTest( // {{{1
		string name_, string pdb_path_, Loop loop_, Size num_solutions_) {

	name = name_;
	pose_from_file(pose, "protocols/kinematic_closure/inputs/" + pdb_path_, core::import_pose::PDB_file);
	original_pose = pose;
	loop = loop_;
	if (loop.cut() == 0) loop.set_cut(loop.midpoint());
	num_solutions = num_solutions_;
	jacobians.resize(num_solutions);
	pivot_picker = PivotPickerOP( new LoopPivots );
	precision = 1e-3;
}

void ClosureTest::perturb_subset( // {{{1
			Pose const &,
			IndexList const &,
			ClosureProblemOP problem) {

	perturb_for_test(problem);
}

void ClosureTest::test_closure(SolutionList solutions) const { // {{{1

	TS_ASSERT_EQUALS(num_solutions, solutions.size());

	foreach (ClosureSolutionCOP solution, solutions) {
		Pose solution_pose(pose);
		solution->apply(solution_pose);

		test_closure(solution_pose);
	}
}

void ClosureTest::test_closure(Pose const & solution_pose) const { // {{{1
	Size left_pivot = index_from_residue(loop.start());
	Size right_pivot = index_from_residue(loop.stop());
	Size start_atom = left_pivot - 3;
	Size stop_atom = right_pivot + 3;

	// Make sure that none of the atoms outside the loop have moved.  In the
	// case that no fold tree is being used, a bug in the closure algorithm
	// will cause atoms outside the designated window to move.

	for (Size index = start_atom; index <= stop_atom; index++) {
		AtomID id = id_from_index(index);
		PointPosition expected = original_pose.xyz(id);
		PointPosition observed = solution_pose.xyz(id);

		if (index <= left_pivot || index >= right_pivot) {
			TS_ASSERT_VECTOR_EQUALS(expected, observed, precision);
		}
	}

	// Make sure that all of the bond lengths are in a reasonable range.  In
	// the case that a fold tree is being used, a bug in the closure algorithm
	// will distort the distance between the two atoms framing the cutpoint.

	AtomID ids[2];
	PointPosition xyzs[2];

	for (Size index = start_atom; index <= start_atom - 1; index++) {
		ids[0] = id_from_index(index);
		ids[1] = id_from_index(index + 1);
		xyzs[0] = solution_pose.xyz(ids[0]);
		xyzs[1] = solution_pose.xyz(ids[1]);

		Real bond_length = xyzs[0].distance(xyzs[1]);

		TR << "index = " << index << endl;
		TR << "ids[0] = " << ids[0] << endl;
		TR << "ids[1] = " << ids[1] << endl;
		TR << "xyzs[0] = " << xyzs[0] << endl;
		TR << "xyzs[1] = " << xyzs[1] << endl;
		TR << "bond_length = " << bond_length << endl;
		TR << endl;

		TS_ASSERT_LESS_THAN(bond_length, 1.6);
		TS_ASSERT_LESS_THAN(1.2, bond_length);
	}
}

void ClosureTest::test_restore( // {{{1
		ClosureProblemCOP problem, SolutionList solutions) const {

	foreach (ClosureSolutionCOP solution, solutions) {
		Pose solution_pose(pose);
		solution->apply(solution_pose);

		// Make sure that the original pose can be restored.

		problem->restore(solution_pose);

		Size num_atoms = 3 * solution_pose.total_residue();

		for (Size index = 1; index <= num_atoms; index++) {
			AtomID id = id_from_index(index);
			PointPosition expected = original_pose.xyz(id);
			PointPosition observed = solution_pose.xyz(id);

			TS_ASSERT_VECTOR_EQUALS(expected, observed, precision);
		}
	}
}

void ClosureTest::test_jacobian(SolutionList solutions) const { // {{{1
	// Make sure the expected jacobians match the observed ones.
	foreach (ClosureSolutionCOP solution, solutions) {
		Real expected_jacobian = jacobians[solution->get_index()];
		Real observed_jacobian = solution->get_jacobian();
		TS_ASSERT_DELTA(expected_jacobian, observed_jacobian, precision);
	}
}

void ClosureTest::test_pickers( // {{{1
		ClosureProblemCOP, SolutionList) const {

	/*
	// Test to make sure the solution pickers work as expected.  The fact that
	// both of the pickers are stochastic is accounted for by making 1000 picks.
	// Within a loose range, the expected distributions should be seen.

	ClosureSolutionCOP solution;

	Size iterations = 1000;
	Real precision = 1.5 * sqrt(iterations);
	Real total_jacobian = 0;

	boost::unordered_map<Size, Size> random_picks;
	boost::unordered_map<Size, Size> balanced_picks;

	foreach (ClosureSolutionCOP solution, solutions) {
		Size index = solution->get_index();
		total_jacobian += solution->get_jacobian();

		random_picks[index] = 0;
		balanced_picks[index] = 0;
	}

	using protocols::kinematic_closure::pick_random_solution;
	using protocols::kinematic_closure::pick_balanced_solution;

	// Make 1000 picks with each picker.  Record the results in hash tables.

	for (int i = 1; i <= iterations; i++) {
		solution = pick_random_solution(solutions);
		random_picks[solution->get_index()] += 1;

		solution = pick_balanced_solution(solutions, solutions);
		balanced_picks[solution->get_index()] += 1;
	}

	// Check to see that the expected distributions were observed.

	foreach (ClosureSolutionCOP solution, solutions) {
		Size index = solution->get_index();

		Size expected_random = iterations / solutions.size();
		Size observed_random = random_picks[index];

		TS_ASSERT_DELTA(expected_random, observed_random, precision);

		Real jacobian_weight = jacobians[index] / total_jacobian;
		Size expected_balanced = iterations * jacobian_weight;
		Size observed_balanced = balanced_picks[index];

		TS_ASSERT_DELTA(expected_balanced, observed_balanced, precision);
	}
	*/
}
// }}}1

AtomID ClosureTest::id_from_index(Size index) const { // {{{1
	return AtomID(((index - 1) % 3) + 1, ((index - 1) / 3) + 1);
}

Size ClosureTest::index_from_residue(Size residue) const { // {{{1
	return 3 * residue - 1;
}
// }}}1

class DebuggingHelper : public ClosureTest { // {{{1

public:

	DebuggingHelper() :
			ClosureTest("debugging_helper", "1ubq.pdb", Loop(2, 6, 3), 4) {

		jacobians[1] = 0.0795;
		jacobians[2] = 0.0879;
		jacobians[3] = 0.1291;
		jacobians[4] = 0.1363;
	}
};

class BigLoopTest : public ClosureTest { // {{{1

public:

	BigLoopTest() :
			ClosureTest("big_loop", "3cla.pdb", Loop(170, 181, 175), 6) {

		jacobians[1] = 0.0032;
		jacobians[2] = 0.0025;
		jacobians[3] = 0.0041;
		jacobians[4] = 0.0029;
		jacobians[5] = 0.0061;
		jacobians[6] = 0.0027;
	}
};

class NoSolutionsTest : public ClosureTest { // {{{1

public:

	NoSolutionsTest() :
			ClosureTest("no_solutions", "1srp.pdb", Loop(313, 319, 316), 0) {
	}

	void perturb_for_test(ClosureProblemOP problem) {
		problem->perturb_phi(313,  -62.040, DEGREES);
		problem->perturb_psi(313,  -31.591, DEGREES);
		problem->perturb_phi(314, -111.245, DEGREES);
  	problem->perturb_psi(314,  155.459, DEGREES);
		problem->perturb_phi(315,   85.594, DEGREES);
  	problem->perturb_psi(315,    2.433, DEGREES);
		problem->perturb_phi(316,  -77.068, DEGREES);
  	problem->perturb_psi(316,   43.512, DEGREES);
		problem->perturb_phi(317,  -90.960, DEGREES);
  	problem->perturb_psi(317,  -66.109, DEGREES);
		problem->perturb_phi(318,  -64.065, DEGREES);
  	problem->perturb_psi(318,  114.326, DEGREES);
		problem->perturb_phi(319,  -91.959, DEGREES);
  	problem->perturb_psi(319,  111.620, DEGREES);
	}
};

class PoseWithLigandTest : public ClosureTest { // {{{1

// If the fold tree is poorly designed, it can get confused by the ligand.

public:

	PoseWithLigandTest() :
			ClosureTest("pose_with_ligand", "1exm.pdb", Loop(289, 300), 2) {

		jacobians[1] = 0.0168;
		jacobians[2] = 0.0172;
	}
};
// }}}1

#endif
