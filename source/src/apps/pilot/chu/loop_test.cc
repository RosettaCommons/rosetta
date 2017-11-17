// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Chu Wang

// Unit Headers
#include <protocols/loops/loops_main.hh>

// Rosetta Headers
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/Pose.hh>
#include <basic/options/option.hh>

#include <protocols/viewer/viewers.hh>

// C++ Headers
#include <iostream>

// option key includes

#include <basic/options/keys/loops.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


void*
my_main( void* )
{

	using namespace basic::options;

	std::cout << "start loops modeling\n";

	core::pose::Pose pose;
	core::import_pose::pose_from_file( pose, option[ OptionKeys::loops::template_pdb ]().name() ); //"input/test_in.pdb", core::import_pose::PDB_file);

	protocols::loops::loops_main( pose );
	return 0;
}

int
main( int argc, char * argv [] )
{
	try {
		// initialize option and random number system
		devel::init( argc, argv );

		protocols::viewer::viewer_main( my_main );
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
