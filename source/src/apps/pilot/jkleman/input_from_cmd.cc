// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    mp_parameters.cc
/// @brief   Checking options for MPframework
/// @details Checking options from database and user input for MPframework
///			     Last Modified: 3/24/14
/// @author  JKLeman (julia.koehler1982@gmail.com)

// App headers
#include <devel/init.hh>

// Project Headers
#include <protocols/membrane/MembraneProteinFactory.hh>
#include <core/conformation/membrane/Exceptions.hh>
#include <core/conformation/membrane/definitions.hh>

// Package Headers
#include <core/types.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <basic/Tracer.hh>

#include <protocols/viewer/viewers.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>

#include <numeric/xyzVector.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>

// C++ Headers
#include <cstdlib>
#include <string>

using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "apps.pilot.jkleman.input_from_cmd" );

using namespace core;

// global data
std::string basename;

///////////////////////////////////////////////////////////////
void create_MPpose(){
	
	using namespace protocols::membrane;
	
	MembraneProteinFactoryOP mpf = new MembraneProteinFactory( false, true );
	core::pose::PoseOP membrane_pose = mpf->create_membrane_pose();
	membrane_pose->dump_pdb("output_pose.pdb");
	
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

