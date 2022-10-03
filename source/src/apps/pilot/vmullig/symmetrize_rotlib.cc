// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/vmullig/symmetrize_rotlib.cc
/// @brief An app that takes a rotamer library that should be symmetric (e.g. a peptoid rotamer library) and generates a symmetrized version.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// devel headers
#include <devel/init.hh>

// protocol headers
#include <protocols/jd2/JobDistributor.hh>
// #include <DUMMYNAMESPACE/DUMMYMOVER.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pack/rotamers/SingleNCAARotamerLibraryCreator.hh>

// utility headers
#include <utility/excn/Exceptions.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <utility/options/OptionCollection.hh>
#include <basic/options/option_macros.hh>

static basic::Tracer TR("symmetrize_rotlib");


OPT_KEY( String, symm_residue_name )

static const std::string errmsg( "Error in symmetrize_rotlib application:  " );

/// @brief Register the relevant options for the --help flag.
void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	NEW_OPT( symm_residue_name, "The residue whose rotamer library to symmetrize.  Required option.", "" );
	// NEW_OPT( output_rotlib, "The output file to generate. Required option.", "" );
}

/// @brief Read from the global options system to get the input and output files.
/// @details Throws error if these options were not provided.
void
get_symm_residue_name(
	std::string & srn
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	runtime_assert_string_msg( option[ symm_residue_name ].user() && !(option[symm_residue_name]().empty()), errmsg + "The -symm_residue_name flag is required!" );
	srn = option[symm_residue_name]();
	// output_file = option[output_rotlib]();
}

/// @brief Entry point for program execution.
int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		register_options();
		devel::init( argc, argv );

		std::string residue_name;
		get_symm_residue_name( residue_name );

		core::chemical::ResidueTypeSetCOP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

		core::chemical::ResidueType const & base_type = rts->name_map( residue_name );

		auto c = core::pack::rotamers::SingleNCAARotamerLibraryCreator();
		c.write_out_ = true;
		c.create(base_type);

		// get_input_and_output_files( infile, outfile );
		// symmetrize_rotlib( infile, outfile );

	} catch ( utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	TR << "Successfully comleted symmetrize_rotlib.  Exiting with code 0 (no errors)." << std::endl;
	return 0;
}
