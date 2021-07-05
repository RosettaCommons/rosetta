// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author jk
#include <string>
#include <sstream>

#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <basic/options/option_macros.hh>

// Utility Headers
#include <utility/vector1.hh>

//Auto Headers

//Protocol Headers
#include <protocols/pockets/PocketGrid.hh>
//#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>

using namespace core;
using namespace basic::options;
using namespace std;
using namespace core::scoring;
using namespace basic::options::OptionKeys;
using namespace conformation;


//OPT_KEY( String, pocket1_fname )
OPT_KEY( Boolean, minipocket )

static basic::Tracer TR( "apps.pilot.david_pocket_compare.main" );

/// General testing code
int main( int argc, char * argv [] ) {

	try {

		//NEW_OPT( pocket1_fname, "pocket", "fname" );
		NEW_OPT( minipocket, "Print smaller set of pocket points as carbons (for ROCS)", false );

		//initializes Rosetta functions
		devel::init(argc, argv);

		for ( core::Size f=1; f <= basic::options::start_files().size(); f++ ) {
			std::string const fname1 = basic::options::start_files().at(f);
			bool const minpock = option [ minipocket ];

			protocols::pockets::TargetPocketGrid pocket1( fname1 );
			std::stringstream out_fname;

			if ( minpock ) {
				out_fname << fname1 << ".min.pdb";
			} else {
				out_fname << fname1 << ".pdb";
			}
			pocket1.dumpTargetPocketsToPDB ( out_fname.str(), minpock);
		}

		TR << "Done!" << std::endl;
	}
catch (utility::excn::Exception const & e ) {
	e.display();
	return -1;
}

	return 0;
}

