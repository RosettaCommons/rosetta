// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    read_pose_from_seq.cc
///
/// @brief   Use the membrane framework reading in a pose from seuqnece. Involves a modified
///    mover that might handle chains differently?
///    Last Modified: 3/27/14
///
/// @author  Rebecca Faye Alford (rfalford12@gmail.com)

// App headers
#include <devel/init.hh>

// Project Headers
#include <protocols/moves/Mover.hh>
#include <protocols/membrane/MembraneUnitTestMover.hh>

// Package Headers
#include <core/types.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <basic/Tracer.hh>

#include <protocols/viewer/viewers.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

using basic::Error;
using basic::Warning;

static basic::Tracer TR( "apps.pilot.ralford.read_mp_pose_from_seq" );

void*
my_main( void* )
{
	using namespace protocols::membrane;
	using namespace protocols::jd2;

	/// how this should work
	/// CreateMPMover( fasta, )
	///

	// Setup MP Unit Testing Mover
	MembraneUnitTestMoverOP mp = new MembraneUnitTestMover();
	JobDistributor::get_instance()->go(mp);

	// Score a Membrane Pose and
	mp->score_lowres();
	mp->score_highres();

	return NULL;
}

///////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

		devel::init(argc, argv);
		protocols::viewer::viewer_main( my_main );

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}

