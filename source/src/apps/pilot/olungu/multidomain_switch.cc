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
#include <devel/domain_assembly/domain_assembly.hh>
#include <basic/options/option.hh>
#include <devel/init.hh>

#include <sstream>
#include <iomanip>
#include <utility/excn/Exceptions.hh>


// option key includes

#include <basic/options/keys/DomainAssembly.OptionKeys.gen.hh>


using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;

//@brief helper function to get viewers to work correctly
void*
my_main( void* ) {
	if ( option[ DomainAssembly::da_setup ] () ) {
		devel::domain_assembly::assemble_domains_setup();
	} else {
		devel::domain_assembly::assemble_domains_optimize();
	}
	return 0;
}


//@brief main function
int
main( int argc, char* argv[] )
{
	try{
	// options, random initialization
	devel::init( argc, argv );

	protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}


