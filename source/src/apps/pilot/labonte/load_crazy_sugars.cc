// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   debug_labontes_current_work.cc
/// @brief  This script tries to load a sugar-containing Pose from a .pdb.  (Yeah, that's it.)
/// @author Labonte <JWLabonte@jhu.edu>


// Package headers
#include <devel/init.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>


using namespace std;
using namespace core;
using namespace pose;
using namespace import_pose;


string const PATH = "/home/labonte/Desktop/labonte/Workspace/For_Frank/";


int
main( int argc, char *argv[] )
{
	try {
		// Initialize core.
		devel::init( argc, argv );

		// Declare variables.
		Pose crazy_sugar;

		// Try to load and output the Pose.
		pose_from_pdb( crazy_sugar, PATH + "fix5-trimerized_chainA_w_LINKs.pdb" );
		cout << crazy_sugar << endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		cerr << "Caught exception: " << e.msg() << endl;
		return -1;
	}
	return 0;
}
