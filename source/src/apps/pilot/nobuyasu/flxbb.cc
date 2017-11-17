// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Nobuyasu Koga

// Unit headers

#include <protocols/flxbb/FlxbbDesign_main.hh>
#include <protocols/viewer/viewers.hh>


// Project headers
//#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/ScoreFunctionFactory.hh>
//#include <core/scoring/constraints/util.hh>
//#include <protocols/jobdist/standard_mains.hh>
//#include <protocols/moves/Mover.hh>


// Utility Headers
#include <basic/options/option.hh>
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>
using basic::Error;
using basic::Warning;

static basic::Tracer TR( "protocols.flxbb" );

// option key includes
//#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/flxbb.OptionKeys.gen.hh>

void*
my_main( void *)
{

	protocols::flxbb::FlxbbDesign_main();
	return 0 ;

}


int
main( int argc, char * argv [] )
{
	try{
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// setup
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		devel::init(argc, argv);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// end of setup
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/*

		core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
		core::scoring::ScoreFunctionOP scorefxn_design = scorefxn;
		core::scoring::ScoreFunctionOP scorefxn_relax = scorefxn;

		if( option[ in::file::fullatom ]() ) {
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *scorefxn_relax );
		} else{
		core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *scorefxn_relax );
		}

		protocols::moves::MoverOP protocol;
		protocol = new protocols::flxbb::FlxbbDesign( scorefxn_design, scorefxn_relax );
		*/

		bool const view( option[ OptionKeys::flxbb::view ] );

		if ( view ) {
			std::cout << "Start viewer mode " << std::endl;
			protocols::viewer::viewer_main( my_main );

		} else {

			protocols::flxbb::FlxbbDesign_main();
		}
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 1;
}


