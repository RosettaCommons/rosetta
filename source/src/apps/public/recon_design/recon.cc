// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/public/recon_design/recon.cc
/// @brief The application file for running RECON multistate design
/// Takes in an XML similar to RosettaScripts and runs RECON from that
/// @author Alex Sevy (alex.sevy@gmail.com)

// Project Headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/VectorPoseJobDistributor.hh>

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

// Register the RECON specific movers that derive from VectorPoseMover/Filter
#include <apps/public/recon_design/init.VectorPoseMoverCreators.ihh>
#include <apps/public/recon_design/init.VectorPoseMoverRegistrators.ihh>
#include <apps/public/recon_design/init.VectorPoseFilterCreators.ihh>
#include <apps/public/recon_design/init.VectorPoseFilterRegistrators.ihh>


// Tracer
static basic::Tracer TR( "apps.public.recon_design.recon" );



/// @brief Runs RECON multistate design as specified in an XML file. Works
/// exactly the same as RosettaScripts, but also allow for use of
/// VectorPoseMovers/Filters in addition to regular movers
int
main( int argc, char * argv [] )
{
	try {
		// not sure what this line does - try it without
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

			protocols::moves::MoverOP mover;//note that this is not instantiated and will crash if the job distributor actually tries to use it.

			if ( !option[ jd2::ntrials ].user() ) {
				// when using rosetta_scripts we want ntrials to be set to 1 if the user forgot to specify. o/w infinite loops might
				// occur.
				option[ jd2::ntrials ].value( 1 );
			}

			// launch the VectorPoseJobDistributor
			protocols::jd2::VectorPoseJobDistributor jd;
			jd.go( mover );
		}

	} catch ( utility::excn::Exception& excn ) {
		excn.display();
		std::exit( 1 );
	}
}

