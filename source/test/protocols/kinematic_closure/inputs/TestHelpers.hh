// -*- 
// mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t 
// -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Header and source file for the TestMover class.
/// @author Kale Kundert (kale.kundert@ucsf.edu)

#ifndef INCLUDED_test_protocols_kinematic_closure_TestHelpers_HH
#define INCLUDED_test_protocols_kinematic_closure_TestHelpers_HH

///////////////////////////////////////////////////////////////////////////////
// Right now, the conservative solution tests are commented out because the 
// algorithm used by the conservative picker counted on the problem making a 
// copy of the input pose.  This was prohibitively slow, so a new algorithm 
// will be needed.  The tests will be reinstated once that algorithm has been 
// written.
///////////////////////////////////////////////////////////////////////////////

#define private public
#define protected public

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/utilities.hh>
#include <protocols/kinematic_closure/ClosureProblem.hh>
#include <protocols/kinematic_closure/ClosureSolution.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/AtomID.hh>

// Protocol headers
#include <protocols/moves/Mover.hh>
#include <protocols/loops/Loop.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <numeric/constants.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <boost/foreach.hpp>
#include <boost/unordered_map.hpp>

// C++ headers
#include <algorithm>
#include <iostream>
#include <sstream>

// Global names {{{1
using namespace std;
using namespace core;

using core::pose::Pose;
using core::pose::PoseOP;
using core::pose::PoseCOP;
using core::chemical::ChemicalManager;
using core::chemical::ResidueTypeSetCAP;
using core::conformation::Conformation;
using core::PointPosition;
using core::id::AtomID;

using protocols::kinematic_closure::ClosureProblemCOP;
using protocols::kinematic_closure::ClosureSolutionCOP;
using protocols::kinematic_closure::ClosureSolution;
using protocols::kinematic_closure::SolutionList;
using protocols::kinematic_closure::IndexList;
using protocols::kinematic_closure::ParameterList;
using protocols::kinematic_closure::setup_fold_tree;

using protocols::moves::Mover;
using protocols::loops::Loop;

using utility::vector1;
using utility::pointer::ReferenceCount;

using numeric::constants::r::pi;

#define foreach BOOST_FOREACH

// Forward declarations {{{1
class TestHelper;
typedef utility::pointer::owning_ptr<TestHelper> TestHelperOP;
typedef utility::pointer::owning_ptr<TestHelper const> TestHelperCOP;

class TestMover;
typedef utility::pointer::owning_ptr<TestMover> TestMoverOP;
typedef utility::pointer::owning_ptr<TestMover const> TestMoverCOP;
// }}}1

// Utility functions
void pose_from_sequence(string const sequence, Pose &pose) { // {{{1
	pose::make_pose_from_sequence(pose, sequence, "fa_standard", false);
}

PoseOP pose_from_sequence(string const sequence) { // {{{1
	PoseOP pose = new Pose();
	pose_from_sequence(sequence, *pose);
	return pose;
}
// }}}1
void write_pdb(Pose const &pose, string path) { // {{{1
	pose.dump_pdb("../../test/protocols/kinematic_closure/" + path);
}

void write_pdb(PoseCOP pose, string path) { // {{{1
	write_pdb(*pose, path);
}

void write_pdbs(Pose const &pose, SolutionList solutions, string prefix) { // {{{1
	string directory = "../../test/protocols/kinematic_closure/";

	foreach (ClosureSolutionCOP solution, solutions) {
		Pose solution_pose(pose);
		solution->apply(solution_pose);

		stringstream path;
		path << directory << prefix << "." << solution->get_index() << ".pdb";

		solution_pose.dump_pdb(path.str());
	}
}
// }}}1

// Utility assertions
void TS_ASSERT_INDICES( // {{{1
		IndexList &indices, Size x, Size y, Size z) {

	TS_ASSERT_EQUALS (indices[1], x);
	TS_ASSERT_EQUALS (indices[2], y);
	TS_ASSERT_EQUALS (indices[3], z);
}

void TS_ASSERT_VECTOR_EQUALS( // {{{1
		PointPosition expected,
		PointPosition observed,
		Real precision) {

	Real distance = expected.distance(observed);
	TS_ASSERT_LESS_THAN(distance, precision);
}

void TS_ASSERT_VECTOR_DIFFERS( // {{{1
		PointPosition expected,
		PointPosition observed,
		Real precision) {

	Real distance = expected.distance(observed);
	TS_ASSERT_LESS_THAN(precision, distance);
}
// }}}1

// Helper classes
class TestHelper : public ReferenceCount { // {{{1

public:
	TestHelper(string, Size, Size, Size, Size, Size);
	string get_name() const { return name; }

public:
	void setup();
	void test_mover(Pose const & pose) const;
	void test_closure(SolutionList solutions) const;
	void test_restore(Pose const & pose) const;
	void test_jacobian(SolutionList solutions) const;
	void test_pickers(ClosureProblemCOP problem, SolutionList solutions) const;

public:
	AtomID id_from_index(Size index) const;
	Size index_from_pivot(Size residue) const;
	bool is_pivot_torsion(Size index) const;

public:
	Pose pose;
	Pose idealized_pose;
	Size num_residues;
	Size num_solutions;
	Loop loop;
	vector1<Size> pivot_residues;
	vector1<PointPosition> atom_xyzs;
	Size conservative_solution;
	vector<Real> jacobians;
	Real precision;
	string name;
};

class TestMover : public Mover { // {{{1

public:
	TestMover(TestHelperOP);
	string get_name() const;

public:
	void apply(Pose &pose);

private:
	TestHelperOP test;
};
// }}}1

TestHelper::TestHelper( // {{{1
		string id,
		Size length,
		Size start,
		Size middle,
		Size end,
		Size solutions) {
	
	name = id;
	num_residues = length;
	num_solutions = solutions;

	loop = Loop(start, end, middle);

	pivot_residues.resize(3);
	pivot_residues[1] = start;
	pivot_residues[2] = middle;
	pivot_residues[3] = end;

	precision = 1e-3;
	atom_xyzs.resize(3 * num_residues);
	jacobians.resize(num_solutions);
}

void TestHelper::setup() { // {{{1
	using core::conformation::Conformation;
	using core::id::AtomID;

	// Initialize the pose with ideal geometry.
	string sequence(num_residues, 'A');
	pose_from_sequence(sequence, pose);
	setup_fold_tree(pose, loop);

	// Create a pose based on the coordinates given in atom_xyzs.
	Size num_atoms = 3 * pose.total_residue();
	Size start_pivot = index_from_pivot(pivot_residues[1]);
	Size end_pivot = index_from_pivot(pivot_residues[3]);

	for (Size index = 1; index <= num_atoms; index++) {
		AtomID id = id_from_index(index);
		pose.set_xyz(id, atom_xyzs[index]);
	}

	// Create another pose where the internal coordinates within the loop have 
	// been idealized.  This is the pose the solution will be applied to, to 
	// guarantee that all the internal degrees of freedom are being set.

	idealized_pose = pose;
	Conformation & conformation = idealized_pose.conformation();
	AtomID ids[4];

	for (Size index = start_pivot; index <= end_pivot - 1; index++) {
		ids[0] = id_from_index(index + 0);
		ids[1] = id_from_index(index + 1);
		conformation.set_bond_length(ids[0], ids[1], 1.4);
	}

	for (Size index = start_pivot; index <= end_pivot - 2; index++) {
		ids[0] = id_from_index(index + 0);
		ids[1] = id_from_index(index + 1);
		ids[2] = id_from_index(index + 2);
		conformation.set_bond_angle(ids[0], ids[1], ids[2], 109.5);
	}

	for (Size index = start_pivot; index <= end_pivot - 3; index++) {
		ids[0] = id_from_index(index + 0);
		ids[1] = id_from_index(index + 1);
		ids[2] = id_from_index(index + 2);
		ids[3] = id_from_index(index + 3);
		conformation.set_torsion_angle(ids[0], ids[1], ids[2], ids[3], 180.0);
	}
}

void TestHelper::test_mover(Pose const & pose) const { // {{{1
	Size num_atoms = 3 * pose.total_residue();
	Size start_pivot = index_from_pivot(pivot_residues[1]);
	Size stop_pivot = index_from_pivot(pivot_residues[3]);

	// Test to make sure the mover actually did something (i.e. broke closure).  
	// This is more of a meta-test than a test of the kinematic closure code.

	for (Size index = 1; index <= num_atoms; index++) {
		AtomID id = id_from_index(index);
		PointPosition open = pose.xyz(id);
		PointPosition closed = atom_xyzs[index];

		if (index > start_pivot && index < stop_pivot) {
			TS_ASSERT_VECTOR_DIFFERS(closed, open, precision);
		}
	}
}

void TestHelper::test_closure(SolutionList solutions) const { // {{{1
	TS_ASSERT_EQUALS(num_solutions, solutions.size());

	foreach (ClosureSolutionCOP solution, solutions) {
		Pose solution_pose(idealized_pose);
		solution->apply(solution_pose);

		Size num_atoms = 3 * solution_pose.total_residue();
		Size start_pivot = index_from_pivot(pivot_residues[1]);
		Size end_pivot = index_from_pivot(pivot_residues[3]);

		// Make sure that none of the atoms outside the loop have moved.  This is 
		// not the strongest possible test, because it's possible that closure 
		// could be broken even if this condition is maintained (think really long 
		// bond lengths).  But it should catch most problems.

		for (Size index = 1; index <= num_atoms; index++) {
			AtomID id = id_from_index(index);
			PointPosition observed = solution_pose.xyz(id);
			PointPosition expected = atom_xyzs[index];

			if (index <= start_pivot || index >= end_pivot) {
				TS_ASSERT_VECTOR_EQUALS(expected, observed, precision);
			}
		}
	}
}

void TestHelper::test_restore(Pose const & pose) const { // {{{1
	Size num_atoms = 3 * pose.total_residue();

	// Make sure that every atom in the pose is exactly where it started.
	for (Size index = 1; index <= num_atoms; index++) {
		AtomID id = id_from_index(index);
		PointPosition observed = pose.xyz(id);
		PointPosition expected = atom_xyzs[index];
		TS_ASSERT_VECTOR_EQUALS(expected, observed, precision);
	}
}

void TestHelper::test_jacobian(SolutionList solutions) const { // {{{1
	// Make sure the expected jacobians match the observed ones.
	foreach (ClosureSolutionCOP solution, solutions) {
		Real expected_jacobian = jacobians[solution->get_index()];
		Real observed_jacobian = solution->get_jacobian();
		TS_ASSERT_DELTA(expected_jacobian, observed_jacobian, precision);
	}
}

void TestHelper::test_pickers( // {{{1
		ClosureProblemCOP problem,
		SolutionList solutions) const {

	// Test to make sure the solution pickers work as expected.  The fact that 
	// two of the three pickers are stochastic is accounted for by making 1000 
	// picks.  Within a loose range, the expected distributions should be seen. 
		
	ClosureSolutionCOP solution;

	Size iterations = 1000;
	Real precision = 1.5 * sqrt(iterations);
	Real total_jacobian = 0;

	boost::unordered_map<Size, Size> random_picks;
	boost::unordered_map<Size, Size> conservative_picks;
	boost::unordered_map<Size, Size> balanced_picks;

	foreach (ClosureSolutionCOP solution, solutions) {
		Size index = solution->get_index();
		total_jacobian += solution->get_jacobian();

		random_picks[index] = 0;
		conservative_picks[index] = 0;
		balanced_picks[index] = 0;
	}

	using protocols::kinematic_closure::pick_random_solution;
	using protocols::kinematic_closure::pick_conservative_solution;
	using protocols::kinematic_closure::pick_balanced_solution;

	// Make 1000 picks with each picker.  Record the results in hash tables.
	
	for (int i = 1; i <= iterations; i++) {
		solution = pick_random_solution(solutions);
		random_picks[solution->get_index()] += 1;

		//solution = pick_conservative_solution(pose, solutions);
		//conservative_picks[solution->get_index()] += 1;

		solution = pick_balanced_solution(solutions, solutions);
		balanced_picks[solution->get_index()] += 1;
	}

	// Check to see that the expected distributions were observed.

	foreach (ClosureSolutionCOP solution, solutions) {
		Size index = solution->get_index();

		Size expected_random = iterations / solutions.size();
		Size observed_random = random_picks[index];

		TS_ASSERT_DELTA(expected_random, observed_random, precision);

		//bool is_conservative_solution = (index == conservative_solution);
		//Size expected_conservative = (is_conservative_solution) ? iterations : 0;
		//Size observed_conservative = conservative_picks[index];

		//TS_ASSERT_EQUALS(expected_conservative, observed_conservative);

		Real jacobian_weight = jacobians[index] / total_jacobian;
		Size expected_balanced = iterations * jacobian_weight;
		Size observed_balanced = balanced_picks[index];

		TS_ASSERT_DELTA(expected_balanced, observed_balanced, precision);
	}
}
// }}}1

AtomID TestHelper::id_from_index(Size index) const { // {{{1
	return AtomID(((index - 1) % 3) + 1, ((index - 1) / 3) + 1);
}

Size TestHelper::index_from_pivot(Size residue) const { // {{{1
	return 3 * residue - 1;
}

bool TestHelper::is_pivot_torsion(Size index) const { // {{{1
	Size start_pivot = index_from_pivot(pivot_residues[1]);
	Size middle_pivot = index_from_pivot(pivot_residues[2]);
	Size end_pivot = index_from_pivot(pivot_residues[3]);

	return (index == start_pivot - 2)  || (index == start_pivot - 1)  ||
	       (index == middle_pivot - 2) || (index == middle_pivot - 1) ||
	       (index == end_pivot - 2)    || (index == end_pivot - 1);
}
// }}}1

TestMover::TestMover(TestHelperOP helper) : test(helper) {} // {{{1

string TestMover::get_name() const { return test->name; } // {{{1

void TestMover::apply(Pose &pose) { // {{{1
	Conformation &conformation = pose.conformation();

	Size num_atoms = 3 * conformation.size();
	Size start_pivot = test->index_from_pivot(test->pivot_residues[1]);
	Size end_pivot = test->index_from_pivot(test->pivot_residues[3]);

	AtomID ids[4];

	// Put all the atoms in the right positions.
	for (Size index = 1; index <= num_atoms; index++) {
		ids[0] = test->id_from_index(index);
		conformation.set_xyz(ids[0], test->atom_xyzs[index]);
	}

	// Change the pivot torsions to create an open pose.
	for (Size index = start_pivot - 2; index <= end_pivot - 1; index++) {
		ids[0] = test->id_from_index(index + 0);
		ids[1] = test->id_from_index(index + 1);
		ids[2] = test->id_from_index(index + 2);
		ids[3] = test->id_from_index(index + 3);

		if (test->is_pivot_torsion(index))
			conformation.set_torsion_angle(ids[0], ids[1], ids[2], ids[3], pi);
	}
}
// }}}1

#endif

