// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/public/apps/public/trRosetta/trRosetta.cc
/// @brief Runs the trRosetta protocol, as described in Yang et al. (2020) Proc. Natl. Acad. Sci. USA
/// 117(3):1496-503.  This is the Python code ported to C++, basically.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// devel headers
#include <devel/init.hh>

// protocol headers
#include <protocols/trRosetta_protocols/movers/trRosettaProtocolMover.hh>


// utility headers
#include <utility/excn/Exceptions.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/tensorflow_manager/util.hh>

#ifdef USE_TENSORFLOW
#include <protocols/jd2/JobDistributor.hh>
#include <basic/options/option.hh>
#endif

static basic::Tracer TR("apps.public.trRosetta.trRosetta");


/// @brief Indicate which commandline flags are relevant to this application.
void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	protocols::trRosetta_protocols::movers::trRosettaProtocolMover::register_options();
}

/// @brief Program entry point.
int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		register_options();
		devel::init( argc, argv );

		TR << "Starting trRosetta application." << std::endl;
		TR << "Application created 21 January 2021 by Vikram K. Mulligan, Flatiron Institute (vmulligan@flatironinstitute.org)." << std::endl;
		TR << "Original trRosetta neural network developed by Jianyi Yang, Ivan Anishchenko, Hahnbeom Park, Zhenling Peng, Sergey Ovchinnikov, and David Baker. (https://doi.org/10.1073/pnas.1914677117)." << std::endl;

#ifndef USE_TENSORFLOW
		utility_exit_with_message( "Could not run trRosetta application!  This version of Rosetta was not compiled "
			"with Tensorflow support.\n\n"
			+ basic::tensorflow_manager::get_tensorflow_compilation_instructions( "trRosetta application" )
		);
#else
		protocols::trRosetta_protocols::movers::trRosettaProtocolMoverOP mover_protocol(
			utility::pointer::make_shared< protocols::trRosetta_protocols::movers::trRosettaProtocolMover >( option )
		);
		protocols::jd2::JobDistributor::get_instance()->go( mover_protocol );
#endif

	} catch ( utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}
