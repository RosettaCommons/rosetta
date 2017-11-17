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

#include <core/types.hh>
#include <devel/init.hh>

#include <core/scoring/methods/SuckerEnergy.hh>
#include <core/pose/Pose.hh>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>


using namespace core;
using namespace scoring;
using namespace ;

void
test_sucker_energy( std::string fname ) {
	methods::SuckerEnergy se;

}


int
main (int argc, char *argv[])
{

	try {


  using namespace basic::options;
  using namespace basic::options::OptionKeys;

	devel::init( argc, argv );

	vector1<file::FileName> files( option[ in::file::s ]() );
	for( size_t i = 1; i <= files.size(); ++i ) {
  	test_sucker_energy( files[i] );
	}

	return 0;


	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
