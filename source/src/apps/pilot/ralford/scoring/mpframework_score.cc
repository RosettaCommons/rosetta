// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    read_mp_pose.cc
///
/// @brief   MP Framework: Read Membrane Pose and Dump Result to PDB
/// @details ast Modified: 3/26/14
///			 
///
/// @author  Rebecca Faye Alford (rfalford12@gmail.com)

// App headers
#include <devel/init.hh>

// Project Headers
#include <protocols/moves/Mover.hh>
#include <protocols/membrane/CreateMembranePoseMover.hh>

#include <core/scoring/ScoreFunction.hh> 
#include <core/scoring/ScoreTypes.hh> 

// Package Headers
#include <core/types.hh>
#include <core/pose/Pose.hh> 

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

static basic::Tracer TR( "apps.pilot.membrane.read_mp_pose" );

void*
my_main( void* )
{
	using namespace protocols::membrane;
	using namespace protocols::jd2;
	
	// Setup MP Unit Testing Mover
	CreateMembranePoseMoverOP mp = new CreateMembranePoseMover();
	JobDistributor::get_instance()->go(mp);
	
	// Score a Membrane Pose and
	core::pose::PoseOP pose = mp->get_membrane_pose();
	
	// Create Sfxn
	core::scoring::ScoreFunctionOP sfxn = new ScoreFunciton();
	sfxn->set_weight( MPPair, 1.0 );
	core::Real total = sfxn->score( pose );
	
	return NULL;
}

///////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		
		devel::init(argc, argv);
		protocols::viewer::viewer_main( my_main );
		
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
	
}

