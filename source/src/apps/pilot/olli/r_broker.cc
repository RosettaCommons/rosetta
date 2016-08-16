// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/olli/r_broker.cc
/// @author Oliver Lange

// keep these headers first for compilation with Visual Studio C++
#include <protocols/jobdist/Jobs.hh>

// Project Headers
#include <protocols/abinitio/BrokerMain.hh>
#include <protocols/checkpoint/Checkpoint.hh>
#include <protocols/viewer/viewers.hh>

// Utility headers
#include <devel/init.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <utility/exit.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <iostream>

#include <utility/vector1.hh>


using namespace protocols;
using namespace abinitio;

void run() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( !option[ in::path::database ].user() ) {
		option[ in::path::database ].def( "/work/olange/minirosetta_database");
	}

	if ( option[ OptionKeys::run::checkpoint ] || option[ OptionKeys::run::checkpoint_interval ].user() ) {
		protocols::checkpoint::checkpoint_with_interval( option[ OptionKeys::run::checkpoint_interval ] );
	}

	Broker_main();
}

void* rBroker_main_local( void* ) {
	try{
		run();
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception : " << std::endl;
		excn.show( std::cerr );
	}
	return 0;
}

int main( int argc, char * argv [] ) {
	try{
		register_options_broker();
		devel::init( argc, argv ); //d
		protocols::viewer::viewer_main( rBroker_main_local );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
