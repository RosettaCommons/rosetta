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


// libRosetta headers
#include <core/types.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/viewer/viewers.hh>
#include <core/pose/Pose.hh>
#include <devel/init.hh>

#include <protocols/recces/RECCES_Mover.hh>
#include <protocols/recces/options/RECCES_Options.hh>
#include <protocols/recces/params/RECCES_Parameters.hh>
#include <protocols/recces/util.hh>
#include <protocols/recces/setup_util.hh>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/recces.OptionKeys.gen.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>

using namespace core::pose;
using namespace basic::options;

using namespace core;
using namespace protocols;
using namespace protocols::moves;
using namespace basic::options::OptionKeys;
using utility::vector1;

///////////////////////////////////////////////////////////////////////////////
void
thermal_sampler()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;
	using namespace protocols::recces;
	using namespace protocols::recces::options;
	using namespace protocols::recces::params;

	RECCES_OptionsOP options( new RECCES_Options );
	options->set_histogram_max( 100.05 ); // sets default
	options->initialize_from_command_line();
	options->set_thermal_sampler_mode( true );
	options->set_blank_score_terms( true );
	options->set_skip_last_accept( true );
	options->set_suppress_sampler_display( true );
	options->set_prefix_each_output_pdb( true );
	options->set_show_more_pose_scores( true );
	options->set_save_scores( true );
	runtime_assert( !options->legacy_turner_mode() );

	PoseOP pose_op( recces_pose_setup( *options ) );
	Pose & pose = *pose_op;
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 600, 600 );
	RECCES_ParametersCOP recces_parameters( new RECCES_Parameters( pose ) );

	RECCES_Mover recces_mover( options );
	recces_mover.set_parameters( recces_parameters );
	recces_mover.set_scorefxn( ( option[ score::weights ].user() ) ? get_score_function() : ScoreFunctionFactory::create_score_function( RNA_HIRES_WTS ) );
	recces_mover.apply( pose );
}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	thermal_sampler();
	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;

		std::cout << std::endl << "Basic usage:  " << argv[0] << "  -s <pdb file> " << std::endl;
		std::cout << std::endl << " Type -help for full slate of options." << std::endl << std::endl;

		option.add_relevant( score::weights );
		option.add_relevant( in::file::s );
		option.add_relevant( in::file::silent );
		option.add_relevant( in::file::tags );
		option.add_relevant( in::file::fasta );
		option.add_relevant( in::file::input_res );
		option.add_relevant( full_model::cutpoint_open );
		option.add_relevant( score::weights );
		option.add_relevant( OptionKeys::recces::out_prefix );
		option.add_relevant( OptionKeys::recces::n_cycle );
		option.add_relevant( OptionKeys::recces::dump_pdb );
		option.add_relevant( OptionKeys::recces::dump_silent );
		option.add_relevant( OptionKeys::recces::out_torsions );
		option.add_relevant( OptionKeys::recces::temps );
		option.add_relevant( OptionKeys::recces::st_weights );
		option.add_relevant( OptionKeys::recces::dump_freq );
		option.add_relevant( OptionKeys::recces::thermal_sampling::sample_residues );
		option.add_relevant( OptionKeys::recces::thermal_sampling::free_residues );
		option.add_relevant( OptionKeys::recces::thermal_sampling::loop_residues );
		option.add_relevant( OptionKeys::recces::thermal_sampling::angle_range_bb );
		option.add_relevant( OptionKeys::recces::thermal_sampling::angle_range_free_bb );
		option.add_relevant( OptionKeys::recces::thermal_sampling::angle_range_chi );
		option.add_relevant( OptionKeys::recces::thermal_sampling::angle_range_free_chi );
		option.add_relevant( OptionKeys::recces::thermal_sampling::chi_stdev );
		option.add_relevant( OptionKeys::recces::thermal_sampling::bb_stdev );
		option.add_relevant( OptionKeys::recces::thermal_sampling::standard_bb_stdev );
		option.add_relevant( OptionKeys::recces::thermal_sampling::setup_base_pair_constraints );

		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		devel::init(argc, argv);

		////////////////////////////////////////////////////////////////////////////
		// end of setup
		////////////////////////////////////////////////////////////////////////////
		protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

