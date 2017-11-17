// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    membrane_sfxn.cc
///
/// @brief   MP Framework: Read Membrane Pose and Dump Result to PDB
/// @details Last Modified: 3/26/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

// App headers
#include <devel/init.hh>

// Project Headers
#include <protocols/moves/Mover.hh>
#include <protocols/membrane/CreateMembranePoseMover.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreTypes.hh>

// New Energy Methods
#include <core/scoring/membrane/MPEnvEnergy.hh>
#include <core/scoring/membrane/MPEnergy.hh>
#include <core/scoring/membrane/MPPairEnergy.hh>
#include <core/scoring/membrane/MPNonHelixPenalty.hh>
#include <core/scoring/membrane/MPTerminiPenalty.hh>
#include <core/scoring/membrane/MPTMProjPenalty.hh>
#include <core/scoring/membrane/MPSpanningPenalty.hh>
#include <core/scoring/membrane/MPSecStructPenalty.hh>

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


	// to be written... would instantiate energy methods and check results?

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

