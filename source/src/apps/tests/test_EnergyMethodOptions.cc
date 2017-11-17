// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/apps/tests/test_EnergyMethodOptions.cc
/// @brief Simply call EnergyMethodOptions::show()
/// @author Andy Watkins (andy.watkins2@gmail.com)

// devel headers
#include <devel/init.hh>

// core headers
#include <core/scoring/methods/EnergyMethodOptions.hh>

// utility headers
#include <utility/excn/Exceptions.hh>

// basic headers
#include <basic/Tracer.hh>

// c++ headers
#include <iostream>


int
main( int argc, char * argv [] )
{
	try {

		devel::init( argc, argv );
		core::scoring::methods::EnergyMethodOptions emo;
		emo.show( std::cout );

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
