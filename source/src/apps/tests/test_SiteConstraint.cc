// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /apps/tests/test_SiteConstraint.cc
/// @brief Application file for testing SiteConstraint.cc
/// @author Morgan Nance <morganlnance@gmail.com>


// Package headers
#include <devel/init.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/constraint_movers/ConstraintSetMover.hh>

#include <utility/excn/Exceptions.hh>

using namespace std;
using namespace core;
using namespace pose;
using namespace import_pose;

string const PATH = "input/";


int
main( int argc, char *argv[] )
{
	try {
		// Initialize core.
		devel::init( argc, argv );

		// Declare variables.
		Pose in_pose1;
		Pose in_pose2;
		protocols::constraint_movers::ConstraintSetMover constraint_setter;

		// Try to load and output the Pose.
		// 3AY4 has three protein chains: A, B, C
		// each with conjugated glycan chains
		// The glycans on chain A and B are the identical and have one branch point
		// One of the glycans on chain C has two branches
		pose_from_file( in_pose1, PATH + "3AY4_modified.pdb" , core::import_pose::PDB_file);
		cout << in_pose1 << endl;
		// Apply the constraints
		constraint_setter.constraint_file( PATH + "3AY4_modified.cst" );
		constraint_setter.apply( in_pose1 );

		// Try to load and output the Pose.
		// 4BJ0 is a carbohydrate-binding module bound to a branched oligosaccharide
		pose_from_file( in_pose2, PATH + "4BJ0.pdb" , core::import_pose::PDB_file);
		cout << in_pose2 << endl;
		// Apply the constraints
		constraint_setter.constraint_file( PATH + "4BJ0.cst" );
		constraint_setter.apply( in_pose2 );

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}
