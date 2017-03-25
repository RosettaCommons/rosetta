// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_test_protocols_kinematic_closure_utilities_HH
#define INCLUDED_test_protocols_kinematic_closure_utilities_HH

// Headers
#include <cxxtest/TestSuite.h>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/import_pose/import_pose.hh>

// Protocol headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/ClosureSolution.hh>

// Utility headers
#include <numeric/constants.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <boost/unordered_map.hpp>

// C++ headers
#include <algorithm>
#include <iostream>
#include <sstream>

// Global names

using std::string;
using std::stringstream;

using core::pose::Pose;
using core::pose::PoseOP;
using core::pose::PoseCOP;

using protocols::kinematic_closure::ClosureSolutionCOP;
using protocols::kinematic_closure::ClosureSolution;
using protocols::kinematic_closure::SolutionList;
using protocols::kinematic_closure::IndexList;
using protocols::kinematic_closure::ParameterList;


// Utility functions
void pose_from_sequence(string const sequence, Pose &pose) { // {{{1
	core::pose::make_pose_from_sequence(pose, sequence, "fa_standard", false);
}

PoseOP pose_from_sequence(string const sequence) { // {{{1
	PoseOP pose( new Pose() );
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

	for ( ClosureSolutionCOP solution: solutions ) {
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
	IndexList &indices, core::Size x, core::Size y, core::Size z) {

	TS_ASSERT_EQUALS (indices[1], x);
	TS_ASSERT_EQUALS (indices[2], y);
	TS_ASSERT_EQUALS (indices[3], z);
}

void TS_ASSERT_VECTOR_EQUALS( // {{{1
	core::PointPosition expected,
	core::PointPosition observed,
	core::Real precision) {

	core::Real distance_between_vectors = expected.distance(observed);
	TS_ASSERT_LESS_THAN(distance_between_vectors, precision);
}

void TS_ASSERT_VECTOR_DIFFERS( // {{{1
	core::PointPosition expected,
	core::PointPosition observed,
	core::Real precision) {

	core::Real distance_between_vectors = expected.distance(observed);
	TS_ASSERT_LESS_THAN(precision, distance_between_vectors);
}
// }}}1

#endif
