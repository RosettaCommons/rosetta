// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/public/stepwise/resample_full_model.cc
/// @brief stepwise-based seq_rebuild (for erraser)
/// @author Caleb Geniesse, geniesse@stanford.edu
/// @author Andy Watkins, amw579@stanford.edu


// protocols
#include <protocols/viewer/viewers.hh>
#include <protocols/stepwise/monte_carlo/util.hh>
#include <protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMasterMover.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/util.hh>


// core
#include <core/types.hh>
#include <devel/init.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/FullModelPoseBuilder.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>

// basic
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>

// utility
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <iostream>
#include <string>

using namespace protocols;
using namespace basic::options;
using namespace basic::options::OptionKeys;





static basic::Tracer TR( "apps.pilot.calebgeniesse.resample_full_model" );

///////////////////////////////////////////////////////////////////////////////
void
resample_full_model_test()
{

	using namespace core::scoring;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::io::silent;
	using namespace core::import_pose;
	using namespace core::import_pose::pose_stream;
	using namespace core::pose;
	using namespace core::pose::rna;
	using namespace core::pose::full_model_info;
	using namespace protocols::stepwise::monte_carlo;

	// setup residue types
	ResidueTypeSetCOP rsd_set;
	rsd_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	// setup score function
	core::scoring::ScoreFunctionOP scorefxn;
	if ( option[ score::weights ].user () ) {
		scorefxn = get_score_function();
	} else {
		scorefxn = ScoreFunctionFactory::create_score_function( RNA_HIRES_WTS );
	}

	// read starting pose
	PDBPoseInputStream input( option[ in::file::s ]() );

	// setup poses
	PoseOP start_pose( new Pose );
	Pose seq_rebuild_pose;


	// iterate over input stream
	input.fill_pose( *start_pose, *rsd_set );
	utility::vector1< PoseOP > input_poses;
	input_poses.push_back( start_pose );
	FullModelPoseBuilder builder;
	builder.set_options( option );
	builder.set_input_poses( input_poses );
	builder.initialize_further_from_options();
	builder.build(); // this will update input_poses[1]
	//fill_full_model_info_from_command_line( start_pose/*, other_ops */);
	( *scorefxn )( *start_pose );

	seq_rebuild_pose = *start_pose;
	protocols::viewer::add_conformation_viewer(seq_rebuild_pose.conformation(), "current", 500, 500);


	// initialize options
	options::StepWiseMonteCarloOptionsOP options( new options::StepWiseMonteCarloOptions );
	options->initialize_from_command_line();
	options->set_erraser( true );
	options->set_enumerate( true );
	options->set_skip_deletions( true );

	// If sample_res is one long, we are in single res mode.
	auto sample_res = const_full_model_info( seq_rebuild_pose ).full_model_parameters()->conventional_to_full( option[ full_model::sample_res ].resnum_and_chain() );
	options->set_output_minimized_pose_list( sample_res.size() == 1 );

	// run StepWiseMasterMover::resample_full_model
	mover::StepWiseMasterMover master_mover( scorefxn, options );
	master_mover.resample_full_model( *start_pose, seq_rebuild_pose, true /*checkpointing_breadcrumbs*/, sample_res );

	// score seq_rebuild pose
	(*scorefxn)(seq_rebuild_pose);

	// show scores of start_pose and full_model_pose
	//if ( option[ show_scores ]() ) {
	std::cout << "\n Score before seq_rebuild:" << std::endl;
	scorefxn->show( std::cout, *start_pose );
	std::cout << "\n Score after seq_rebuild:" << std::endl;
	scorefxn->show( std::cout, seq_rebuild_pose );
	//}

	// dump pdb -- should figure out file name based on inputs
	seq_rebuild_pose.dump_pdb( "SEQ_REBUILD.pdb" );


	// write out to silent file
	//std::string tag = "S_0";
	//stepwise::monte_carlo::output_to_silent_file( tag, silent_file_out, full_model_pose );


}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	resample_full_model_test();
	exit( 0 );
}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

		// help
		std::cout << std::endl << "Basic usage:  " << argv[0] << "  -in::file::silent <silent file>" << std::endl;
		std::cout << std::endl << " Type -help for full slate of options." << std::endl << std::endl;

		// options

		// setup
		devel::init(argc, argv);
		option[ OptionKeys::chemical::patch_selectors ].push_back( "VIRTUAL_BASE" );
		option[ OptionKeys::chemical::patch_selectors ].push_back( "VIRTUAL_RIBOSE" );
		option[ OptionKeys::chemical::patch_selectors ].push_back( "VIRTUAL_RNA_RESIDUE" );
		option[ OptionKeys::chemical::patch_selectors ].push_back( "TERMINAL_PHOSPHATE" );

		// run
		protocols::viewer::viewer_main( my_main );
		exit( 0 );

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
