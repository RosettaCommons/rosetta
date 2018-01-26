// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/jadolfbr/rosetta_scripts_jd3.cc
/// @brief First implementation of a JD3 app for RosettaScripts.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifdef USEMPI
#include <mpi.h>
#endif

// devel headers
#include <devel/init.hh>

// protocol headers
#include <protocols/jd3/JobDistributor.hh>
#include <protocols/jd3/JobDistributorFactory.hh>
#include <protocols/jd3/util.hh>

#include <protocols/rosetta_scripts/RosettaScriptsJobQueen.hh>
#include <protocols/rosetta_scripts/util.hh>

// utility headers
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>


static basic::Tracer TR("rosetta_scripts_jd3");

using namespace protocols::rosetta_scripts;
using namespace protocols::jd3;

int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		devel::init( argc, argv );

		//Enables quick one of running as same behavior as rosetta_scripts application.

		if ( option[ parser::info ].user() ) { // If the -parser::info option is used, just print information about the requested component(s) and exit.
			protocols::rosetta_scripts::print_information( option[ parser::info ]() );
		} else if ( option[ parser::output_schema ].user() ) {
			protocols::rosetta_scripts::save_schema( option[ parser::output_schema ] );
		} else if ( ! option[ in::file::job_definition_file ].user() && ! option[ parser::protocol ].user() ) { // Just print a template script and exit if no input script is provided.

			protocols::rosetta_scripts::print_template_script();
			protocols::jd3::print_job_template();
			TR << "\n\nExample of Basic <Options>: \n" <<
				"<parser__protocol value=\"my_rosetta_script.xml\"/>" << std::endl <<
				"<parser__script_vars value=\"branch=1 cartmin=0 rounds=15 min_rings=0\" />\n\n" << std::endl << std::endl;

		} else { // If an input script has been provided, then we're not printing a template script and exiting.

			//View does not make sense for JD3 - I think more would need to be done to get that compatible.  I think watching a single process should be possible,
			// but not sure what it would take to get this to work.  For now, we skip it.

			protocols::jd3::JobDistributorOP jd = protocols::jd3::JobDistributorFactory::create_job_distributor();

			protocols::jd3::JobQueenOP queen( new protocols::rosetta_scripts::RosettaScriptsJobQueen );
			jd->go( queen );
		}

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
