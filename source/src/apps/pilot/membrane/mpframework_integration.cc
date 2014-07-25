// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		src/apps/pilot/membrane/mpframework_integration.cc
///
/// @brief 		Membrane Framework Integration
/// @details    Relax with membrane highres sfxn modifications and over
///				detailed embeddings
///
/// @author 	Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified: 6/15/14

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <protocols/membrane/AddMembraneMover.hh> 

// Project Headers
#include <protocols/jd2/util.hh>
#include <protocols/jd2/JobDistributor.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <iostream>

using namespace protocols;

/// @brief Main method
int
main( int argc, char * argv [] )
{
	try {

		using namespace protocols::membrane;

		protocols::jd2::register_options();
	
		devel::init(argc, argv);

		AddMembraneMoverOP add_mem = new AddMembraneMover(); 
		protocols::jd2::JobDistributor::get_instance()->go( add_mem );

		return 0; 

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

