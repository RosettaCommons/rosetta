// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    visualize_membrane.cc
///
/// @brief   Read in a pose into the membranew framework and
///          dump out the pose with the lipid bilayer present (as VRT atoms)
/// @details last Modified: 4/4/14
/// @author  JKLeman (julia.koehler1982@gmail.com)

// App headers
#include <devel/init.hh>

// Project Headers
#include <protocols/moves/Mover.hh>
#include <protocols/membrane/CreateMembranePoseMover.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>


// Package Headers
#include <apps/benchmark/performance/init_util.hh>
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

static thread_local basic::Tracer TR( "apps.pilot.membrane.visualize_membrane" );

void create_MPpose(){
 
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::membrane;

	// TODO: in::file::membrane_input needs to be taken out once bogus embedding stuff is implemented
	// and once the framework can deal with missing files (spanfile and embedding stuff)
	option[ OptionKeys::in::membrane ]( true );
	option[ OptionKeys::in::file::membrane_input ]( "input_spanfiles.txt" );
	if ( ! option[ out::membrane_pdb_thickness ].user()){
			option[ OptionKeys::out::membrane_pdb_thickness ]( 30 );
	}
	
	// create MP
	MembraneProteinFactoryOP mpf = new MembraneProteinFactory( false, true );
	core::pose::PoseOP membrane_pose = mpf->create_membrane_pose();
	membrane_pose->dump_pdb("membrane_pose.pdb");
    
}


///////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
        
		devel::init(argc, argv);
        
		// create MP pose
		create_MPpose();
        
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
    
}
