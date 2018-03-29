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
		Pose in_pose;

		// Try to load and output the Pose.
		pose_from_file( in_pose, PATH + "protein_with_glycans.pdb" , core::import_pose::PDB_file);
		cout << in_pose << endl;

		protocols::constraint_movers::ConstraintSetMover constraint_setter;
		constraint_setter.constraint_file( PATH + "glycan_to_protein_and_glycan.cst" );
		constraint_setter.apply( in_pose );

	} catch (utility::excn::Exception const & e ) {
		cerr << "Caught exception: " << e.msg() << endl;
		return -1;
	}
	return 0;
}
