// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license.
// (c) The Rosetta software is developed by the contributing members of the
// (c) Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org.
// (c) Questions about this can be addressed to University of Washington UW
// (c) TechTransfer, email: license@u.washington.edu.

/// @file   test_non-canonical_connectivities.cc
/// @brief  This is an integration test app for testing general non-canonical connectivities..
/// @author Labonte <JWLabonte@jhu.edu>

// includes
#include <devel/init.hh>

//#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>

#include <utility/excn/Exceptions.hh>

#include <iostream>


using namespace std;
using namespace core;
using namespace import_pose;
using namespace pose;


string const INPATH = "input/";
string const OUTPATH = "output/";


int
main( int argc, char *argv[] )
{
	using namespace std;

	try {
		// Initialize core.
		devel::init( argc, argv );

		// Import test poses.
		Pose lactam;
		pose_from_file( lactam, INPATH + "lactam.pdb" , core::import_pose::PDB_file);

		cout << "Intra-peptide lactam PDB file" << endl;
		cout << lactam << endl;
		cout << "5: " << lactam.residue( 5 ).name() << endl;
		cout << "11: " << lactam.residue( 11 ).name() << endl;

		lactam.dump_pdb( OUTPATH + "lactam.pdb" );
	} catch ( utility::excn::EXCN_Base const & e ) {
		cerr << "caught exception " << e.msg() << endl;
		return -1;
	}
	return 0;
}
