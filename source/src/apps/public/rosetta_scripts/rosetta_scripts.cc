// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/public/rosetta_scripts/rosetta_scripts.cc
/// @brief The application file for rosetta_scripts, aka jd2_scripting or the parser
/// @author Sarel Fleishman (sarelf@u.washington.edu)

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
#include <protocols/moves/Mover.fwd.hh>

#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>

#include <utility/vector1.hh>
#include <utility/excn/EXCN_Base.hh>
#include <basic/Tracer.hh>

// C++ headers

void* my_main( void *);

void*
my_main( void *)
{
	protocols::moves::MoverOP mover;//note that this is not instantiated and will crash if the job distributor actually tries to use it. That means that this can only be used with parser=true
	protocols::jd2::JobDistributor::get_instance()->go(mover);
	return 0 ;
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

		bool const view( option[ parser::view ] );
		protocols::moves::MoverOP mover;//note that this is not instantiated and will crash if the job distributor actually tries to use it.
																		//That means that this can only be used with parser=true
		option[ jd2::dd_parser ].value( true ); // So here we fix that. jd2_parser app makes no sense without this option=true
		if( !option[ jd2::ntrials ].user() )
	// when using rosteta_scripts we want ntrials to be set to 1 if the user forgot to specify. o/w infinite loops might
	// occur. We don't want ntrials to be set as default to 1, b/c other protocols might want it to behave differently
			option[ jd2::ntrials ].value( 1 );

		if( view )
			protocols::viewer::viewer_main( my_main );
		else{
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
	} catch( utility::excn::EXCN_Base& excn ) {
		basic::Error()
			<< "ERROR: Exception caught by rosetta_scripts application:"
			<< excn << std::endl;
		std::exit( 1 );
	}
}

