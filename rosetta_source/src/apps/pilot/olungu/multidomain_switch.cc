// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief

#include <protocols/viewer/viewers.hh>
#include <devel/domain_assembly/domain_assembly_setup.hh>
// AUTO-REMOVED #include <devel/domain_assembly/domain_assembly.hh>
#include <basic/options/option.hh>
#include <devel/init.hh>

#include <sstream>
#include <iomanip>


// option key includes

#include <basic/options/keys/DomainAssembly.OptionKeys.gen.hh>


using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;

// // protocol specific options
// namespace DomainAssembly{
//   BooleanOptionKey const da_setup( "DomainAssembly::da_setup" );
//   FileOptionKey const da_setup_output_pdb( "DomainAssembly::da_setup_output_pdb" );
// 	FileOptionKey const da_setup_option_file( "DomainAssembly::da_setup_option_file" );
// 	FileOptionKey const start_pdb( "DomainAssembly::da_start_pdb" );
// 	FileOptionKey const linker_file( "DomainAssembly::da_linmer_file" );
// 	IntegerOptionKey const start_pdb_num( "DomainAssembly::da_start_pdb_num" );
// 	IntegerOptionKey const nruns( "DomainAssembly::da_nruns" );
// }


//@brief helper function to get viewers to work correctly
void*
my_main( void* ) {
	if ( option[ DomainAssembly::da_setup ] () ) {
		assemble_domains_setup();
	} else {
		assemble_domains_optimize();
	}
	return 0;
}

//
//@brief main function
int
main( int argc, char* argv[] )
{
	// options, random initialization
	devel::init( argc, argv );

	protocols::viewer::viewer_main( my_main );
	return 0;
}


