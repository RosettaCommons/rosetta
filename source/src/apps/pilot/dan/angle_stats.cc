// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief currently just prints out omega angles annotated for cis for the whole pose
/// @author Daniel J. Mandell

#include <protocols/viewer/viewers.hh>
#include <core/io/pdb/pose_io.hh>

#include <core/pose/Pose.hh>
#include <devel/init.hh>
#include <basic/Tracer.hh>

// option key includes
//#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <basic/options/option.hh>
#include <cmath>

static thread_local basic::Tracer TR( "anglestats" );

////////////////////////////////////////////////////////
void *
anglestats_local( void* ) {
	using namespace core;
	using namespace basic::options;

	pose::Pose pose;
	core::import_pose::pose_from_pdb( pose, option[ OptionKeys::in::file::s ]().vector().front() );
	for (Size i=2; i<= pose.total_residue()-1; i++) {
		if (std::abs(pose.omega(i)) < 150) {
			TR << "cis ";
		}
		TR << "omega for residue " << i << ": " << pose.omega(i) << std::endl;
	}
	return 0;
}

////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
	// options, random initialization
	devel::init( argc, argv );
	protocols::viewer::viewer_main( anglestats_local );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

