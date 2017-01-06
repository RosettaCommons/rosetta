// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Legacy version of recces app, to be merged with thermal_sampler by December 2016

// libRosetta headers
#include <devel/init.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/recces/RECCES_Mover.hh>
#include <protocols/recces/options/RECCES_Options.hh>
#include <protocols/recces/params/RECCES_Parameters.hh>
#include <protocols/recces/setup_util.hh>
#include <protocols/viewer/viewers.hh>


// C++ headers
#include <string>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/recces.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Exception handling
#include <utility/excn/Exceptions.hh>

static THREAD_LOCAL basic::Tracer TR( "recces_turner" );

using namespace core;
using namespace basic::options::OptionKeys;
using namespace basic::options;

//////////////////////////////////////////////////////////////////////////////
void
MC_run() {

	using namespace core::pose;
	using namespace core::scoring;
	using namespace protocols::recces;
	using namespace protocols::recces::params;
	using namespace protocols::recces::options;

		TR << TR.Red << "This app recces_turner will soon be deprecated. Instead use recces executable with flags -seq1, etc." << std::endl;

	RECCES_OptionsOP recces_options( new RECCES_Options );
	recces_options->initialize_from_command_line();
	runtime_assert( recces_options->legacy_turner_mode() ); // must supply -seq1 (maybe -seq2 too, ...)

	Pose pose( *recces_pose_setup( *recces_options ) );
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 600, 600 );
	RECCES_ParametersCOP recces_parameters( new RECCES_Parameters( pose ) );

	RECCES_Mover recces_mover( recces_options );
	recces_mover.set_parameters( recces_parameters );
	recces_mover.set_scorefxn( ( option[ score::weights ].user() ) ? get_score_function() : ScoreFunctionFactory::create_score_function( RNA_HIRES_WTS ) );

	recces_mover.apply( pose );

}
//////////////////////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	MC_run();
	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}
//////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	using namespace core;

	try {
		option.add_relevant( OptionKeys::recces::temps );
		option.add_relevant( OptionKeys::recces::st_weights );
    option.add_relevant( OptionKeys::recces::n_cycle );
		option.add_relevant( OptionKeys::recces::seq1 );
		option.add_relevant( OptionKeys::recces::seq2 );
    option.add_relevant( OptionKeys::recces::dump_pdb );

		devel::init ( argc, argv );
		protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

