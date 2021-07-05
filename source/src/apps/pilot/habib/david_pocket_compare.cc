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
#include <fstream>
#include <string>

#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <basic/options/option_macros.hh>

// Utility Headers

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


OPT_KEY( String, pocket1_fname )
OPT_KEY( String, pocket2_fname )

static basic::Tracer TR( "apps.pilot.david_pocket_compare.main" );

/// General testing code
int main( int argc, char * argv [] ) {

	try{

		NEW_OPT( pocket1_fname, "pocket1_fname", "fname" );
		NEW_OPT( pocket2_fname, "pocket2_fname", "fname" );

		//initializes Rosetta functions
		devel::init(argc, argv);

		std::string const fname1 ( option[ pocket1_fname ] );
		std::string const fname2 ( option[ pocket2_fname ] );

		TR << "Starting pocket compare" << std::endl;

		protocols::pockets::TargetPocketGrid pocket1( fname1 );
		protocols::pockets::TargetPocketGrid pocket2( fname2 );


		// call function to compare template and comparison Pockets, report score
		core::Real d1 = pocket1.get_pocket_distance( pocket2 );
		TR << "Distance is: " << d1 << std::endl;

		//utility::io::ozstream fout;
		//fout.open("distance.txt", std::ios::out);
		//fout << d1 << std::endl;
		//fout.close();
		//  fout.clear();

		TR << "Done!" << std::endl;
	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;

}

