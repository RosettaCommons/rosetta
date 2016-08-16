// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/devel/balanced_kic/assertion_helpers.hh
/// @brief  Provide several convenient assertion functions.
/// @author Kale Kundert
///
/// The CxxTest Suite provides a number of primitive assertion macros, but 
/// these are inconvenient to use for higher-order data structures like vectors 
/// or poses.  This header provides functions (not macros) that compose the 
/// primitive assertions into more convenient units.  Note that the functions 
/// are given macro-like names, for the sake of consistency.

#ifndef INCLUDED_devel_balanced_kic_assertion_helpers_HH
#define INCLUDED_devel_balanced_kic_assertion_helpers_HH

#include <cxxtest/TestSuite.h>
#include <devel/balanced_kic/algorithms.hh>

#include <core/pose/Pose.hh>
#include <core/id/NamedAtomID.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <sstream>

using namespace core;
using namespace devel::balanced_kic;
using namespace devel::balanced_kic::algorithms;

void TS_ASSERT_VECTOR(
		Coordinate &vector, Real x, Real y, Real z, Real epsilon) {

	TS_ASSERT_DELTA (vector[1], x, epsilon);
	TS_ASSERT_DELTA (vector[2], y, epsilon);
	TS_ASSERT_DELTA (vector[3], z, epsilon);
}

void TS_ASSERT_VECTOR(
		pose::Pose &pose, Size index, std::string name, 
		Real x, Real y, Real z, Real epsilon) {

	id::NamedAtomID id (name, index);
	numeric::xyzVector<Real> position = pose.xyz(id);

	TS_ASSERT_DELTA (position.x(), x, epsilon);
	TS_ASSERT_DELTA (position.y(), y, epsilon);
	TS_ASSERT_DELTA (position.z(), z, epsilon);
}

void TS_ASSERT_COPIED_POSE(
		pose::Pose const *first_pose,
		pose::Pose const *second_pose,
		bool expected_to_be_identical=true ) {

	std::stringstream first_pdb, second_pdb;

	first_pose->dump_pdb(first_pdb);
	second_pose->dump_pdb(second_pdb);

	// Ensure that both poses generate the same PDB file, but that they also 
	// have different addresses in memory.

	TS_ASSERT_DIFFERS (first_pose, second_pose);
	if (expected_to_be_identical) {
		TS_ASSERT_EQUALS (first_pdb.str(), second_pdb.str());
	}
}

#endif
