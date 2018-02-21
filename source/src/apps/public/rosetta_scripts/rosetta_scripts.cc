// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/public/rosetta_scripts/rosetta_scripts.cc
/// @brief The application file for rosetta_scripts, aka jd2_scripting or the parser
/// @author Sarel Fleishman (sarelf@u.washington.edu)
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- Added template script generation and "-parser:info" flag (in-app help for RosettaScripts).

// Project Headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/MPIWorkPoolJobDistributor.hh>
#include <protocols/jd2/MPIFileBufJobDistributor.hh>
#include <protocols/jd2/BOINCJobDistributor.hh>
#include <protocols/viewer/viewers.hh>
#include <core/types.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>

#include <devel/init.hh>
#include <basic/options/option.hh>

// Utility Headers

// Unit Headers
#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/FilterFactory.hh>
#include <protocols/moves/MoverFactory.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <core/pack/task/operation/ResLvlTaskOperationFactory.hh>

#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

#include <basic/Tracer.hh>

#include <fstream>

// Tracer
static basic::Tracer TR( "apps.public.rosetta_scripts.rosetta_scripts" );

// FUNCTION PROTOTYPES
void* my_main( void *);



// FUNCTION DECLARATIONS
void*
my_main( void *)
{
	protocols::moves::MoverOP mover;//note that this is not instantiated and will crash if the job distributor actually tries to use it. That means that this can only be used with parser=true
	protocols::jd2::JobDistributor::get_instance()->go(mover);
	return nullptr ;
}



/// @details dock_design_scripting provides an xml-based scripting capability
/// to run rosetta movers and filters defined in a text file provided by the
/// user. A full documentation of dock_design_scripting is available at:
/// manual_doxygen/applications/app_dock_design.dox
int
main( int argc, char * argv [] )
{
	try{
		protocols::abinitio::ClassicAbinitio::register_options();
		// setup random numbers and options
		devel::init(argc, argv);
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		if ( option[ parser::info ].user() ) { // If the -parser::info option is used, just print information about the requested component(s) and exit.
			protocols::rosetta_scripts::print_information( option[ parser::info ]() );
		} else if ( option[ parser::output_schema ].user() ) {
			protocols::rosetta_scripts::save_schema( option[ parser::output_schema ] );
		} else if ( ! option[ parser::protocol ].user() ) { // Just print a template script and exit if no input script is provided.
			protocols::rosetta_scripts::print_template_script();
		} else { // If an input script has been provided, then we're not printing a template script and exiting.
			bool const view( option[ parser::view ] );
			protocols::moves::MoverOP mover;//note that this is not instantiated and will crash if the job distributor actually tries to use it.

			if ( !option[ jd2::ntrials ].user() ) {
				// when using rosetta_scripts we want ntrials to be set to 1 if the user forgot to specify. o/w infinite loops might
				// occur.
				option[ jd2::ntrials ].value( 1 );
			}

			if ( view ) {
				protocols::viewer::viewer_main( my_main );
			} else {
#ifdef BOINC
				protocols::jd2::BOINCJobDistributor::get_instance()->go( mover );
#else
#ifdef USEMPI
				protocols::jd2::MPIFileBufJobDistributor::get_instance()->go( mover );
#else
				protocols::jd2::JobDistributor::get_instance()->go( mover );
#endif
#endif
			}
		}
	} catch (utility::excn::Exception& excn ) {
		basic::Error()
			<< "Exception caught by rosetta_scripts application:"
			<< excn << std::endl;
		std::exit( 1 );
	}
}

