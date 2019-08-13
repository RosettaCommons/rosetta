// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/vmullig/make_mainchain_potential.cc
/// @brief An app to make a mainchain potential, using MM or QM methods, for an NCAA with an already-generated sidechain potential
/// (Dunbrack library).
/// @details I'll convert this to a JD3 app for parallelization later.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// devel headers
#include <devel/init.hh>

// protocol headers
#include <protocols/mainchain_potential/GenerateMainchainPotential.hh>
#include <protocols/mainchain_potential/GenerateMainchainPotentialOptions.hh>

// utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/pointer/memory.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option_macros.hh>

static basic::Tracer TR("make_mainchain_potential");


void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	protocols::mainchain_potential::GenerateMainchainPotentialOptions::register_options();

	//option.add_relevant( in::file::s );
	//option.add_relevant( in::file::l );
}

int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		register_options();
		devel::init( argc, argv );

		protocols::mainchain_potential::GenerateMainchainPotentialOptionsOP options( utility::pointer::make_shared< protocols::mainchain_potential::GenerateMainchainPotentialOptions >( true /*Initialize from global options system.*/ ) );
		protocols::mainchain_potential::GenerateMainchainPotential generator(options);
		generator.run();
		generator.write_last_generated_to_disk();
	} catch ( utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}
