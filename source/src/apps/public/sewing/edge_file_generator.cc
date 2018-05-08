// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/guffysl/edge_file_generator.cc
/// @brief App to test new SewGraphGenerator class
/// @author guffysl (guffy@email.unc.edu)

// devel headers
#include <devel/init.hh>

// protocol headers
#include <protocols/sewing/hashing/EdgeMapGenerator.hh>

// utility headers
#include <utility/excn/Exceptions.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <utility/options/OptionCollection.hh>

static basic::Tracer TR("edge_file_generator");


int
main( int argc, char * argv [] )
{
	try {
		devel::init( argc, argv );
		TR << "devel::init successful" << std::endl;
		TR << "Created SewGraphGenerator" << std::endl;
		protocols::sewing::hashing::EdgeMapGeneratorOP generator = protocols::sewing::hashing::EdgeMapGeneratorOP(new protocols::sewing::hashing::EdgeMapGenerator);
		generator->initialize_from_options();
		generator->generate_edge_file();

	} catch ( utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
