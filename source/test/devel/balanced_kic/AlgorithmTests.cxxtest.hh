// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/devel/balanced_kic/KinematicMover.cxxtest.hh
/// @brief  Test the algorithms that drive the KinematicMover.
/// @author Kale Kundert
///
/// Each of the tests in this suite exercises one of the helper functions used 
/// by the Kinematic Mover.  In most of the tests, the helper function in 
/// question is provided with one set of input data.  The resulting output data 
/// is then exhaustively checked, either against hand-calculated results or 
/// results from Dan Mandell's KinematicMover implementation.  
/// 
/// The input data is taken from the ``loop.pdb`` file located in this 
/// directory.  It is a 5-residue segment of an ideal alpha-helix, build using 
/// Pymol, with pivots in the 2, 3, and 4 positions.  Two closure solutions 
/// exist for this configuration, but one has a near-infinite jacobian weight.  
/// In the future, I'd like to find a different input structure that has two 
/// closure solutions with more reasonable jacobian weights.

#ifndef INCLUDED_devel_balanced_kic_algorithm_tests_CXXTEST_HH
#define INCLUDED_devel_balanced_kic_algorithm_tests_CXXTEST_HH

#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/devel/balanced_kic/assertion_helpers.hh>
#include <test/devel/balanced_kic/dummy_perturber.hh>

#include <devel/balanced_kic/algorithms.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>
#include <boost/foreach.hpp>

using namespace core;
using namespace devel::balanced_kic;
using namespace devel::balanced_kic::algorithms;

class AlgorithmTests : public CxxTest::TestSuite {

public:

	void setUp() { // {{{1
		core_init();

		import_pose::pose_from_file(pose, "devel/balanced_kic/loop.pdb", core::import_pose::PDB_file);

		define_closure_problem(problem, pose, 2, 3, 4);
		solve_closure_problem(problem, solutions);
		apply_all_closure_solutions(solutions);
		normalize_all_closure_solutions(solutions);
	}

	void tearDown() { // {{{1
		problem = ClosureProblem();
		solutions = SolutionList();
	}
	// }}}1

	// Tested algorithms.

	void test_define_closure_problem() { // {{{1

		// Test Loop Definition {{{2

		TS_ASSERT (problem.has_definition);
		TS_ASSERT_COPIED_POSE (&pose, &problem.pose);

		TS_ASSERT_EQUALS (problem.first_residue, 2);
		TS_ASSERT_EQUALS (problem.num_residues, 3);
		TS_ASSERT_EQUALS (problem.pivot_atoms[1], 5);
		TS_ASSERT_EQUALS (problem.pivot_atoms[2], 8);
		TS_ASSERT_EQUALS (problem.pivot_atoms[3], 11);

		// Test Atom Positions {{{2

		TS_ASSERT_EQUALS (problem.atom_xyzs.size(), 15);

		TS_ASSERT_VECTOR (problem.atom_xyzs[ 1],  -2.614, -3.357, -1.789,  1e-3);
		TS_ASSERT_VECTOR (problem.atom_xyzs[ 2],  -1.165, -3.357, -1.789,  1e-3);
		TS_ASSERT_VECTOR (problem.atom_xyzs[ 3],  -0.641, -1.928, -1.789,  1e-3);
                                 
		TS_ASSERT_VECTOR (problem.atom_xyzs[ 4],  -1.201, -1.087, -2.662,  1e-3);
		TS_ASSERT_VECTOR (problem.atom_xyzs[ 5],  -0.790,  0.300, -2.759,  1e-3);
		TS_ASSERT_VECTOR (problem.atom_xyzs[ 6],  -0.945,  0.987, -1.410,  1e-3);
                                 
		TS_ASSERT_VECTOR (problem.atom_xyzs[ 7],  -2.097,  0.789, -0.767,  1e-3);
		TS_ASSERT_VECTOR (problem.atom_xyzs[ 8],  -2.366,  1.389,  0.525,  1e-3);
		TS_ASSERT_VECTOR (problem.atom_xyzs[ 9],  -1.293,  0.981,  1.524,  1e-3);
                                 
		TS_ASSERT_VECTOR (problem.atom_xyzs[10],  -0.975, -0.315,  1.564,  1e-3);
		TS_ASSERT_VECTOR (problem.atom_xyzs[11],   0.031, -0.829,  2.472,  1e-3);
		TS_ASSERT_VECTOR (problem.atom_xyzs[12],   1.355, -0.116,  2.243,  1e-3);
                                 
		TS_ASSERT_VECTOR (problem.atom_xyzs[13],   1.757,  0.008,  0.976,  1e-3);
		TS_ASSERT_VECTOR (problem.atom_xyzs[14],   2.999,  0.667,  0.627,  1e-3);
		TS_ASSERT_VECTOR (problem.atom_xyzs[15],   3.010,  2.087,  1.175,  1e-3);

		// Test Bond Lengths {{{2

		TS_ASSERT_EQUALS (problem.bond_lengths.size(), 15);

		TS_ASSERT_DELTA (problem.bond_lengths[ 1], 1.449, 1e-3);
		TS_ASSERT_DELTA (problem.bond_lengths[ 2], 1.522, 1e-3);
		TS_ASSERT_DELTA (problem.bond_lengths[ 3], 1.335, 1e-3);
		TS_ASSERT_DELTA (problem.bond_lengths[ 4], 1.449, 1e-3);
		TS_ASSERT_DELTA (problem.bond_lengths[ 5], 1.521, 1e-3);
		TS_ASSERT_DELTA (problem.bond_lengths[ 6], 1.334, 1e-3);
		TS_ASSERT_DELTA (problem.bond_lengths[ 7], 1.449, 1e-3);
		TS_ASSERT_DELTA (problem.bond_lengths[ 8], 1.521, 1e-3);
		TS_ASSERT_DELTA (problem.bond_lengths[ 9], 1.335, 1e-3);
		TS_ASSERT_DELTA (problem.bond_lengths[10], 1.449, 1e-3);
		TS_ASSERT_DELTA (problem.bond_lengths[11], 1.521, 1e-3);
		TS_ASSERT_DELTA (problem.bond_lengths[12], 1.335, 1e-3);
		TS_ASSERT_DELTA (problem.bond_lengths[13], 1.448, 1e-3);
		TS_ASSERT_DELTA (problem.bond_lengths[14], 1.522, 1e-3);
		TS_ASSERT_DELTA (problem.bond_lengths[15], 8.370, 1e-3);

		// Test Bond Angles {{{2
		
		TS_ASSERT_EQUALS (problem.bond_angles.size(), 15);

		TS_ASSERT_DELTA (problem.bond_angles[ 1],  47.782, 1e-3);
		TS_ASSERT_DELTA (problem.bond_angles[ 2], 110.137, 1e-3);
		TS_ASSERT_DELTA (problem.bond_angles[ 3], 116.547, 1e-3);
		TS_ASSERT_DELTA (problem.bond_angles[ 4], 121.828, 1e-3);
		TS_ASSERT_DELTA (problem.bond_angles[ 5], 110.102, 1e-3);
		TS_ASSERT_DELTA (problem.bond_angles[ 6], 116.629, 1e-3);
		TS_ASSERT_DELTA (problem.bond_angles[ 7], 121.894, 1e-3);
		TS_ASSERT_DELTA (problem.bond_angles[ 8], 110.076, 1e-3);
		TS_ASSERT_DELTA (problem.bond_angles[ 9], 116.608, 1e-3);
		TS_ASSERT_DELTA (problem.bond_angles[10], 121.895, 1e-3);
		TS_ASSERT_DELTA (problem.bond_angles[11], 110.097, 1e-3);
		TS_ASSERT_DELTA (problem.bond_angles[12], 116.648, 1e-3);
		TS_ASSERT_DELTA (problem.bond_angles[13], 121.941, 1e-3);
		TS_ASSERT_DELTA (problem.bond_angles[14], 110.111, 1e-3);
		TS_ASSERT_DELTA (problem.bond_angles[15],  42.340, 1e-3);

		// Test Dihedral Angles {{{2

		TS_ASSERT_EQUALS (problem.torsion_angles.size(), 15);

		TS_ASSERT_DELTA (problem.torsion_angles[ 1], 331.434, 1e-3);
		TS_ASSERT_DELTA (problem.torsion_angles[ 2], 313.043, 1e-3);
		TS_ASSERT_DELTA (problem.torsion_angles[ 3], 179.965, 1e-3);
		TS_ASSERT_DELTA (problem.torsion_angles[ 4], 302.933, 1e-3);
		TS_ASSERT_DELTA (problem.torsion_angles[ 5], 313.051, 1e-3);
		TS_ASSERT_DELTA (problem.torsion_angles[ 6], 180.029, 1e-3);
		TS_ASSERT_DELTA (problem.torsion_angles[ 7], 302.976, 1e-3);
		TS_ASSERT_DELTA (problem.torsion_angles[ 8], 313.027, 1e-3);
		TS_ASSERT_DELTA (problem.torsion_angles[ 9], 180.001, 1e-3);
		TS_ASSERT_DELTA (problem.torsion_angles[10], 302.957, 1e-3);
		TS_ASSERT_DELTA (problem.torsion_angles[11], 313.015, 1e-3);
		TS_ASSERT_DELTA (problem.torsion_angles[12], 179.979, 1e-3);
		TS_ASSERT_DELTA (problem.torsion_angles[13], 303.066, 1e-3);
		TS_ASSERT_DELTA (problem.torsion_angles[14], 327.347, 1e-3);
		TS_ASSERT_DELTA (problem.torsion_angles[15], 348.880, 1e-3);

		// }}}2
	
	}

	void test_perturb_loop_coordinates() { // {{{1

		DummyPerturber perturber;
		bool exhausted;

		TS_ASSERT_DIFFERS (problem.bond_lengths[1], -1)
		TS_ASSERT_DIFFERS (problem.bond_angles[1], -1)
		TS_ASSERT_DIFFERS (problem.torsion_angles[1], -1)

		perturb_loop_coordinates(problem, perturber, exhausted);

		TS_ASSERT_EQUALS (problem.bond_lengths[1], -1)
		TS_ASSERT_EQUALS (problem.bond_angles[1], -1)
		TS_ASSERT_EQUALS (problem.torsion_angles[1], -1)

		TS_ASSERT (exhausted);

	}

	void test_solve_closure_problem() { // {{{1

		// Test Loop Definition {{{2

		TS_ASSERT_EQUALS (solutions.size(), 2);

		TS_ASSERT (solutions[1].has_internal_coordinates);
		TS_ASSERT_EQUALS (solutions[1].first_residue, 2);
		TS_ASSERT_EQUALS (solutions[1].num_residues, 3);
		TS_ASSERT_EQUALS (solutions[1].pivot_atoms[1], 5);
		TS_ASSERT_EQUALS (solutions[1].pivot_atoms[2], 8);
		TS_ASSERT_EQUALS (solutions[1].pivot_atoms[3], 11);
		TS_ASSERT_COPIED_POSE (&problem.pose, &solutions[1].pose, false);

		TS_ASSERT (solutions[2].has_internal_coordinates);
		TS_ASSERT_EQUALS (solutions[2].first_residue, 2);
		TS_ASSERT_EQUALS (solutions[2].num_residues, 3);
		TS_ASSERT_EQUALS (solutions[2].pivot_atoms[1], 5);
		TS_ASSERT_EQUALS (solutions[2].pivot_atoms[2], 8);
		TS_ASSERT_EQUALS (solutions[2].pivot_atoms[3], 11);
		TS_ASSERT_COPIED_POSE (&problem.pose, &solutions[2].pose, false);

		// Test Bond Lengths {{{2

		TS_ASSERT_EQUALS (solutions[1].bond_lengths.size(), 15);
		TS_ASSERT_EQUALS (solutions[2].bond_lengths.size(), 15);

		TS_ASSERT_DELTA (solutions[1].bond_lengths[ 1], 1.449, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_lengths[ 2], 1.522, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_lengths[ 3], 1.335, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_lengths[ 4], 1.450, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_lengths[ 5], 1.522, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_lengths[ 6], 1.334, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_lengths[ 7], 1.450, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_lengths[ 8], 1.522, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_lengths[ 9], 1.335, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_lengths[10], 1.449, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_lengths[11], 1.521, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_lengths[12], 1.335, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_lengths[13], 1.449, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_lengths[14], 1.522, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_lengths[15], 8.370, 1e-3);

		TS_ASSERT_DELTA (solutions[2].bond_lengths[ 1], 1.449, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_lengths[ 2], 1.522, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_lengths[ 3], 1.335, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_lengths[ 4], 1.450, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_lengths[ 5], 1.522, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_lengths[ 6], 1.334, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_lengths[ 7], 1.450, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_lengths[ 8], 1.522, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_lengths[ 9], 1.335, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_lengths[10], 1.449, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_lengths[11], 1.521, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_lengths[12], 1.335, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_lengths[13], 1.449, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_lengths[14], 1.522, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_lengths[15], 8.370, 1e-3);

		// Test Bond Angles {{{2
		
		TS_ASSERT_EQUALS (solutions[1].bond_angles.size(), 15);
		TS_ASSERT_EQUALS (solutions[2].bond_angles.size(), 15);

		TS_ASSERT_DELTA (solutions[1].bond_angles[ 1],  47.782, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_angles[ 2], 110.137, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_angles[ 3], 116.547, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_angles[ 4], 121.828, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_angles[ 5], 110.102, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_angles[ 6], 116.629, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_angles[ 7], 121.894, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_angles[ 8], 110.076, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_angles[ 9], 116.608, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_angles[10], 121.895, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_angles[11], 110.097, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_angles[12], 116.648, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_angles[13], 121.941, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_angles[14], 110.111, 1e-3);
		TS_ASSERT_DELTA (solutions[1].bond_angles[15],  42.340, 1e-3);

		TS_ASSERT_DELTA (solutions[2].bond_angles[ 1],  47.782, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_angles[ 2], 110.137, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_angles[ 3], 116.547, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_angles[ 4], 121.828, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_angles[ 5], 110.102, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_angles[ 6], 116.629, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_angles[ 7], 121.894, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_angles[ 8], 110.076, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_angles[ 9], 116.608, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_angles[10], 121.895, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_angles[11], 110.097, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_angles[12], 116.648, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_angles[13], 121.941, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_angles[14], 110.111, 1e-3);
		TS_ASSERT_DELTA (solutions[2].bond_angles[15],  42.340, 1e-3);

		// Test Dihedral Angles {{{2

		TS_ASSERT_EQUALS (solutions[1].torsion_angles.size(), 15);
		TS_ASSERT_EQUALS (solutions[2].torsion_angles.size(), 15);

		TS_ASSERT_DELTA (solutions[1].torsion_angles[ 1], 331.434, 1e-3);
		TS_ASSERT_DELTA (solutions[1].torsion_angles[ 2], 313.043, 1e-3);
		TS_ASSERT_DELTA (solutions[1].torsion_angles[ 3], 179.965, 1e-3);
		TS_ASSERT_DELTA (solutions[1].torsion_angles[ 4], 302.933, 1e-3);
		TS_ASSERT_DELTA (solutions[1].torsion_angles[ 5], 313.051, 1e-3);
		TS_ASSERT_DELTA (solutions[1].torsion_angles[ 6], 180.029, 1e-3);
		TS_ASSERT_DELTA (solutions[1].torsion_angles[ 7], 302.976, 1e-3);
		TS_ASSERT_DELTA (solutions[1].torsion_angles[ 8], 313.027, 1e-3);
		TS_ASSERT_DELTA (solutions[1].torsion_angles[ 9], 180.001, 1e-3);
		TS_ASSERT_DELTA (solutions[1].torsion_angles[10], 302.957, 1e-3);
		TS_ASSERT_DELTA (solutions[1].torsion_angles[11], 313.015, 1e-3);
		TS_ASSERT_DELTA (solutions[1].torsion_angles[12], 179.979, 1e-3);
		TS_ASSERT_DELTA (solutions[1].torsion_angles[13], 303.066, 1e-3);
		TS_ASSERT_DELTA (solutions[1].torsion_angles[14], 327.347, 1e-3);
		TS_ASSERT_DELTA (solutions[1].torsion_angles[15], 348.880, 1e-3);

		TS_ASSERT_DELTA (solutions[2].torsion_angles[ 1], 331.434, 1e-3);
		TS_ASSERT_DELTA (solutions[2].torsion_angles[ 2], 313.043, 1e-3);
		TS_ASSERT_DELTA (solutions[2].torsion_angles[ 3], 179.965, 1e-3);
		TS_ASSERT_DELTA (solutions[2].torsion_angles[ 4], 298.814, 1e-3);
		TS_ASSERT_DELTA (solutions[2].torsion_angles[ 5], 324.279, 1e-3);
		TS_ASSERT_DELTA (solutions[2].torsion_angles[ 6], 180.029, 1e-3);
		TS_ASSERT_DELTA (solutions[2].torsion_angles[ 7], 293.379, 1e-3);
		TS_ASSERT_DELTA (solutions[2].torsion_angles[ 8], 321.557, 1e-3);
		TS_ASSERT_DELTA (solutions[2].torsion_angles[ 9], 180.001, 1e-3);
		TS_ASSERT_DELTA (solutions[2].torsion_angles[10], 292.011, 1e-3);
		TS_ASSERT_DELTA (solutions[2].torsion_angles[11], 315.598, 1e-3);
		TS_ASSERT_DELTA (solutions[2].torsion_angles[12], 179.979, 1e-3);
		TS_ASSERT_DELTA (solutions[2].torsion_angles[13], 303.066, 1e-3);
		TS_ASSERT_DELTA (solutions[2].torsion_angles[14], 327.347, 1e-3);
		TS_ASSERT_DELTA (solutions[2].torsion_angles[15], 348.880, 1e-3);

		// }}}2
		
	}

	void test_perform_sanity_check() { // {{{1

		TS_ASSERT_EQUALS (perform_sanity_check(solutions[1]), true);
		solutions[1].torsion_angles[5] = 1e100;
		TS_ASSERT_EQUALS (perform_sanity_check(solutions[1]), false);

	}

	void test_apply_closure_solution() { // {{{1

		// Test Atom Positions {{{2

		CoordinateList &atom_xyzs_1 = solutions[1].atom_xyzs;
		CoordinateList &atom_xyzs_2 = solutions[2].atom_xyzs;

		// First Solution

		TS_ASSERT_EQUALS (atom_xyzs_1.size(), 15);

		TS_ASSERT_VECTOR (atom_xyzs_1[ 1],  -2.614, -3.357, -1.789,  1e-3);
		TS_ASSERT_VECTOR (atom_xyzs_1[ 2],  -1.165, -3.357, -1.789,  1e-3);
		TS_ASSERT_VECTOR (atom_xyzs_1[ 3],  -0.641, -1.928, -1.789,  1e-3);

		TS_ASSERT_VECTOR (atom_xyzs_1[ 4],  -1.201, -1.087, -2.662,  1e-3);
		TS_ASSERT_VECTOR (atom_xyzs_1[ 5],  -0.790,  0.300, -2.759,  1e-3);
		TS_ASSERT_VECTOR (atom_xyzs_1[ 6],  -0.945,  0.987, -1.410,  1e-3);
		
		TS_ASSERT_VECTOR (atom_xyzs_1[ 7],  -2.097,  0.789, -0.767,  1e-3);
		TS_ASSERT_VECTOR (atom_xyzs_1[ 8],  -2.366,  1.389,  0.525,  1e-3);
		TS_ASSERT_VECTOR (atom_xyzs_1[ 9],  -1.293,  0.981,  1.524,  1e-3);

		TS_ASSERT_VECTOR (atom_xyzs_1[10],  -0.975, -0.315,  1.564,  1e-3);
		TS_ASSERT_VECTOR (atom_xyzs_1[11],   0.031, -0.829,  2.472,  1e-3);
		TS_ASSERT_VECTOR (atom_xyzs_1[12],   1.355, -0.116,  2.243,  1e-3);

		TS_ASSERT_VECTOR (atom_xyzs_1[13],   1.757,  0.008,  0.976,  1e-3);
		TS_ASSERT_VECTOR (atom_xyzs_1[14],   2.999,  0.667,  0.627,  1e-3);
		TS_ASSERT_VECTOR (atom_xyzs_1[15],   3.010,  2.087,  1.175,  1e-3);

		// Second Solution

		TS_ASSERT_EQUALS (atom_xyzs_2.size(), 15);

		TS_ASSERT_VECTOR (atom_xyzs_2[ 1],  -2.614, -3.357, -1.789,  1e-3);
		TS_ASSERT_VECTOR (atom_xyzs_2[ 2],  -1.165, -3.357, -1.789,  1e-3);
		TS_ASSERT_VECTOR (atom_xyzs_2[ 3],  -0.641, -1.928, -1.789,  1e-3);

		TS_ASSERT_VECTOR (atom_xyzs_2[ 4],  -1.201, -1.087, -2.662,  1e-3);
		TS_ASSERT_VECTOR (atom_xyzs_2[ 5],  -0.790,  0.300, -2.759,  1e-3);
		TS_ASSERT_VECTOR (atom_xyzs_2[ 6],  -1.040,  1.013, -1.438,  1e-3);

		TS_ASSERT_VECTOR (atom_xyzs_2[ 7],  -2.128,  0.648, -0.758,  1e-3);
		TS_ASSERT_VECTOR (atom_xyzs_2[ 8],  -2.476,  1.252,  0.513,  1e-3);
		TS_ASSERT_VECTOR (atom_xyzs_2[ 9],  -1.441,  0.886,  1.567,  1e-3);

		TS_ASSERT_VECTOR (atom_xyzs_2[10],  -0.959, -0.358,  1.524,  1e-3);
		TS_ASSERT_VECTOR (atom_xyzs_2[11],   0.031, -0.829,  2.472,  1e-3);
		TS_ASSERT_VECTOR (atom_xyzs_2[12],   1.355, -0.116,  2.243,  1e-3);

		TS_ASSERT_VECTOR (atom_xyzs_2[13],   1.757,  0.008,  0.976,  1e-3);
		TS_ASSERT_VECTOR (atom_xyzs_2[14],   2.999,  0.667,  0.627,  1e-3);
		TS_ASSERT_VECTOR (atom_xyzs_2[15],   3.010,  2.087,  1.175,  1e-3);

		// Test Pose Position {{{2

		pose::Pose &pose_1 = solutions[1].pose;
		pose::Pose &pose_2 = solutions[2].pose;

		// First Solution

		TS_ASSERT_VECTOR (pose_1, 1, "N",   -2.614, -3.357, -1.789,  1e-3);
		TS_ASSERT_VECTOR (pose_1, 1, "CA",  -1.165, -3.357, -1.789,  1e-3);
		TS_ASSERT_VECTOR (pose_1, 1, "C",   -0.641, -1.928, -1.789,  1e-3);

		TS_ASSERT_VECTOR (pose_1, 2, "N",   -1.201, -1.087, -2.662,  1e-3);
		TS_ASSERT_VECTOR (pose_1, 2, "CA",  -0.790,  0.300, -2.759,  1e-3);
		TS_ASSERT_VECTOR (pose_1, 2, "C",   -0.945,  0.987, -1.410,  1e-3);
		
		TS_ASSERT_VECTOR (pose_1, 3, "N",   -2.097,  0.789, -0.767,  1e-3);
		TS_ASSERT_VECTOR (pose_1, 3, "CA",  -2.366,  1.389,  0.525,  1e-3);
		TS_ASSERT_VECTOR (pose_1, 3, "C",   -1.293,  0.981,  1.524,  1e-3);

		TS_ASSERT_VECTOR (pose_1, 4, "N",   -0.975, -0.315,  1.564,  1e-3);
		TS_ASSERT_VECTOR (pose_1, 4, "CA",   0.031, -0.829,  2.472,  1e-3);
		TS_ASSERT_VECTOR (pose_1, 4, "C",    1.355, -0.116,  2.243,  1e-3);

		TS_ASSERT_VECTOR (pose_1, 5, "N",    1.757,  0.008,  0.976,  1e-3);
		TS_ASSERT_VECTOR (pose_1, 5, "CA",   2.999,  0.667,  0.627,  1e-3);
		TS_ASSERT_VECTOR (pose_1, 5, "C",    3.010,  2.087,  1.175,  1e-3);

		// Second Solution

		TS_ASSERT_VECTOR (pose_2, 1, "N",   -2.614, -3.357, -1.789,  1e-3);
		TS_ASSERT_VECTOR (pose_2, 1, "CA",  -1.165, -3.357, -1.789,  1e-3);
		TS_ASSERT_VECTOR (pose_2, 1, "C",   -0.641, -1.928, -1.789,  1e-3);

		TS_ASSERT_VECTOR (pose_2, 2, "N",   -1.201, -1.087, -2.662,  1e-3);
		TS_ASSERT_VECTOR (pose_2, 2, "CA",  -0.790,  0.300, -2.759,  1e-3);
		TS_ASSERT_VECTOR (pose_2, 2, "C",   -1.040,  1.013, -1.438,  1e-3);

		TS_ASSERT_VECTOR (pose_2, 3, "N",   -2.128,  0.648, -0.758,  1e-3);
		TS_ASSERT_VECTOR (pose_2, 3, "CA",  -2.476,  1.252,  0.513,  1e-3);
		TS_ASSERT_VECTOR (pose_2, 3, "C",   -1.441,  0.886,  1.567,  1e-3);

		TS_ASSERT_VECTOR (pose_2, 4, "N",   -0.959, -0.358,  1.524,  1e-3);
		TS_ASSERT_VECTOR (pose_2, 4, "CA",   0.031, -0.829,  2.472,  1e-3);
		TS_ASSERT_VECTOR (pose_2, 4, "C",    1.355, -0.116,  2.243,  1e-3);

		TS_ASSERT_VECTOR (pose_2, 5, "N",    1.757,  0.008,  0.976,  1e-3);
		TS_ASSERT_VECTOR (pose_2, 5, "CA",   2.999,  0.667,  0.627,  1e-3);
		TS_ASSERT_VECTOR (pose_2, 5, "C",    3.010,  2.087,  1.175,  1e-3);

		// }}}2

	}

	void test_normalize_closure_solution() { // {{{1

		TS_ASSERT_DELTA (solutions[1].jacobian, 1.550, 1e-3);

		// The jacobian for the second solution is on the order of 10^100.  Since I 
		// don't actually know if this result is correct, I feel uncomfortable 
		// explicitly testing for it.

		//TS_ASSERT_DELTA (solutions[2].jacobian, 1, 1e-3);

	}

	void test_chained_solution_list() { // {{{1

		SolutionList a, b;
		SolutionList const &c = a, &d = b;
		ChainedSolutionList solutions(c, d);

		for (int i = 0; i < 3; i++) {
			ClosureSolution x;
			ClosureSolution y;

			x.id = i;
			y.id = i + 3;

			a.push_back(x);
			b.push_back(y);
		}

		Size expected_id = 0;

		BOOST_FOREACH (ClosureSolution const &solution, solutions) {
			TS_ASSERT_EQUALS (solution.id, expected_id);
			expected_id += 1;
		}

		TS_ASSERT_EQUALS (expected_id, 6);

	}

	void test_pick_random_solution() { // {{{1

		ClosureSolution solution;
		float iterations = 1e4;
		utility::vector1<int> pick_counter(2);

		for (int i = 1; i <= iterations; i++) {
			pick_random_solution(solutions, solution);
			pick_counter[solution.id]++;
		}

		TS_ASSERT_DELTA (pick_counter[1] / iterations, 0.5, 0.1)
		TS_ASSERT_DELTA (pick_counter[2] / iterations, 0.5, 0.1)

	}

	void test_pick_balanced_solution() { // {{{1
 
		ClosureSolution solution;
		Real total_jacobian = solutions[1].jacobian + solutions[2].jacobian;
		float iterations = 1e4;
		utility::vector1<int> pick_counter(2);

		for (int i = 1; i <= iterations; i++) {
			pick_balanced_solution(solutions, solutions, solution);
			pick_counter[solution.id]++;
		}

		TS_ASSERT_DELTA (
				pick_counter[1] / iterations, 
				solutions[1].jacobian / total_jacobian, 0.1)

		TS_ASSERT_DELTA (
				pick_counter[2] / iterations, 
				solutions[2].jacobian / total_jacobian, 0.1)

	}
	// }}}1

	// Untested algorithms.

	void idealize_loop_coordinates() {	// {{{1

		//ClosureProblem problem;
		//idealize_loop_coordinates(problem);

	} // }}}1

private:

	pose::Pose pose;
	ClosureProblem problem;
	SolutionList solutions;

};

#endif


