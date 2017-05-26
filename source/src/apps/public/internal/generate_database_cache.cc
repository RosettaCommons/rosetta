// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief Generate any of the pre-cached files from the database
/// @details Intended to be called in a single-processor context only, with an empty commandline
/// @file Rocco Moretti (rmorettiase@gmail.com)


#include <core/pack/dunbrack/RotamerLibrary.hh>

#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <utility/exit.hh>

#include <core/init/score_function_corrections.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/mistakes.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

static THREAD_LOCAL basic::Tracer TR( "apps.public.generate_database_cache" );

void
generate_dunbrack_binaries() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// I'm using scoping here to make sure the RotamerLibraries are unloaded shortly after they're made

	// Generate the default Rotamer library
	{
		TR.Debug << "Making sure the default Dunbrack binaries are generated" << std::endl;
		core::pack::dunbrack::RotamerLibrary* rotamer_library(  core::pack::dunbrack::RotamerLibrary::get_instance() );
		TR << "Generated default Dunbrack binary cache: " << rotamer_library->get_binary_name() << std::endl;
	}

	// Generate the beta Rotamer Library
	{
		TR.Debug << "Making sure the beta Dunbrack binaries are generated" << std::endl;
		utility::options::OptionCollection beta_options( basic::options::option );
		beta_options[ corrections::beta ].set_value( "true" );
		core::init::init_score_function_corrections( beta_options );

		core::pack::dunbrack::RotamerLibrary beta_rotamer_library( beta_options );
		TR << "Generated Dunbrack binary cache for -beta: " << beta_rotamer_library.get_binary_name() << std::endl;
	}

	// Generate the talaris Rotamer Library
	{
		TR.Debug << "Making sure the talaris Dunbrack binaries are generated" << std::endl;
		utility::options::OptionCollection talaris_options( basic::options::option );
		talaris_options[ corrections::restore_talaris_behavior ].set_value( "true" );
		core::init::init_score_function_corrections( talaris_options );

		core::pack::dunbrack::RotamerLibrary talaris_rotamer_library( talaris_options );
		TR << "Generated Dunbrack binary cache for -restore_talaris_behavior: " << talaris_rotamer_library.get_binary_name() << std::endl;
	}

	// Generate the score12 Rotamer Library
	{
		TR.Debug << "Making sure the score12 Dunbrack binaries are generated." << std::endl;
		utility::options::OptionCollection score12_options( basic::options::option );
		score12_options[ mistakes::restore_pre_talaris_2013_behavior ].set_value( "true" );
		core::init::init_score_function_corrections( score12_options );

		core::pack::dunbrack::RotamerLibrary score12_rotamer_library( score12_options );
		TR << "Generated Dunbrack binary cache for -restore_pre_talaris_2013_behavior: " << score12_rotamer_library.get_binary_name() << std::endl;
	}

	// Generate the score12 -correct Rotamer Library
	{
		TR.Debug << "Making sure the score12 -correct Dunbrack binaries are generated." << std::endl;
		utility::options::OptionCollection score12correct_options( basic::options::option );
		score12correct_options[ mistakes::restore_pre_talaris_2013_behavior ].set_value( "true" );
		score12correct_options[ corrections::correct ].set_value( "true" );
		core::init::init_score_function_corrections( score12correct_options );

		core::pack::dunbrack::RotamerLibrary score12correct_rotamer_library( score12correct_options );
		TR << "Generated Dunbrack binary cache for -restore_pre_talaris_2013_behavior -correct: " << score12correct_rotamer_library.get_binary_name() << std::endl;
	}

}

int
main( int argc, char * argv [] )
{
	try {

		// initialize core
		devel::init(argc, argv);

		generate_dunbrack_binaries();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

