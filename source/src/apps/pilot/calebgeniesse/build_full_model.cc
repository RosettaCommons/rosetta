// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/calebgeniesse/build_full_model.cc
/// @brief Build in missing residues with stepwise::monte_carlo::build_full_model
/// @author Caleb Geniesse, geniesse@stanford.edu


// protocols
#include <protocols/farna/RNA_LoopCloser.hh>
#include <protocols/farna/RNA_LoopCloser.fwd.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/stepwise/monte_carlo/util.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/file_util.hh>
#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>


// core
#include <core/types.hh>
#include <core/init/init.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/VariantType_mappings.cc>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/conformation/Residue.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
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

// utility
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <iostream>
#include <string>

using namespace protocols;
using namespace basic::options;
using namespace basic::options::OptionKeys;

OPT_KEY( Boolean, virtualize_built )
OPT_KEY( Boolean, show_scores )
OPT_KEY( Boolean, dump_pdb )


static THREAD_LOCAL basic::Tracer TR( "apps.pilot.calebgeniesse.build_full_model" );

///////////////////////////////////////////////////////////////////////////////
void
build_full_model_test()
{

	using namespace core::scoring;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::io::silent;
	using namespace core::import_pose::pose_stream;
	using namespace core::pose;
	using namespace core::pose::rna;
	using namespace core::pose::full_model_info;

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

	// setup silent_files
	if ( !option[ in::file::silent ].user() ) {
		utility_exit_with_message( "Error: Must supply a silent file! Use -in::file::silent <silent file>" );
	}
	std::string const silent_file_out = option[ out::file::silent ]();
	if ( option[ out::overwrite ]() ) {
		stepwise::modeler::remove_silent_file_if_it_exists( silent_file_out );
	}

	// setup input stream
	PoseInputStreamOP input;
	if ( option[ in::file::tags ].user() ) {
		input = PoseInputStreamOP( new SilentFilePoseInputStream(
			option[ in::file::silent ](),
			option[ in::file::tags ]()
			));
	} else {
		input = PoseInputStreamOP( new SilentFilePoseInputStream(
			option[ in::file::silent ]()
			));
	}

	// setup poses
	Pose start_pose, full_model_pose;
	utility::vector1< PoseOP > other_ops;

	// setup movers
	farna::RNA_LoopCloserOP rna_loop_closer( farna::RNA_LoopCloserOP( new farna::RNA_LoopCloser() ) );
	rna_loop_closer->set_verbose( true );

	// iterate over input stream
	while ( input->has_another_pose() ) {

		// fill/score start_pose
		input->fill_pose( start_pose, *rsd_set );
		other_ops = nonconst_full_model_info( start_pose ).other_pose_list();
		stepwise::setup::fill_full_model_info_from_command_line( start_pose, other_ops );
		std::string const tag = tag_from_pose( start_pose );
		(*scorefxn)(start_pose);

		// build full_model_pose
		TR << "Building full model for pose: " << TR.Cyan << tag << TR.Reset << std::endl;
		if ( other_ops.size() ) {
			TR << TR.Red << "[WARNING] other_pose_list.size() != 0, full_model_pose may be garbage!!! " << TR.Reset << std::endl;
		}
		stepwise::monte_carlo::build_full_model( start_pose, full_model_pose );

		// compare full_model_info of start pose and full_model_pose
		FullModelInfo const & start_info = const_full_model_info( start_pose );
		FullModelInfo const & full_model_info = const_full_model_info( full_model_pose );

		// virtualize built residues
		for ( Size i = 1; i <= full_model_info.res_list().size(); ++i ) {

			Size const & full_model_res = full_model_info.res_list()[i];
			Size other_idx = start_info.get_idx_for_other_pose_with_residue( full_model_res );

			if ( start_info.res_list().has_value( full_model_res ) || other_idx ) {

				// re-apply variants of residue in start_pose, to residue in full_model_pose
				Size const & sub_res = start_info.full_to_sub( full_model_res );
				Residue const & start_rsd = other_idx ? other_ops[other_idx]->residue( sub_res ) : start_pose.residue( sub_res );
				utility::vector1< std::string > variant_types = start_rsd.type().properties().get_list_of_variants();

				for ( Size j = 1; j <= variant_types.size(); ++j ) {
					VariantType const & variant_type = start_rsd.type().properties().get_variant_from_string( variant_types[j] );
					add_variant_type_to_pose_residue( full_model_pose, variant_type, full_model_res );
				}

			} else {

				// close cutpoints
				other_idx = start_info.get_idx_for_other_pose_with_residue( full_model_res + 1 );
				if ( start_info.res_list().has_value( full_model_res + 1 ) || other_idx ) {
					TR << "Closing loop at residue " << TR.Cyan << full_model_res << TR.Reset << std::endl;
					add_variant_type_to_pose_residue( full_model_pose, CUTPOINT_LOWER, full_model_res );
					add_variant_type_to_pose_residue( full_model_pose, CUTPOINT_UPPER, full_model_res + 1 );
					rna_loop_closer->apply( full_model_pose ); //, full_model_res );
				}

				// apply virtual_rna_residue variant to residue in full_model_pose
				if ( option[ virtualize_built ]() ) {
					TR << "Virtualizing residue " << TR.Cyan << full_model_res << TR.Reset << std::endl;
					apply_virtual_rna_residue_variant_type( full_model_pose, full_model_res, true );
					remove_variant_type_from_pose_residue( full_model_pose, CUTPOINT_LOWER, full_model_res );
					remove_variant_type_from_pose_residue( full_model_pose, CUTPOINT_UPPER, full_model_res + 1 );
				}

			}
		}

		// update pdb info/ score full_model_pose
		update_pdb_info_from_full_model_info( full_model_pose );
		(*scorefxn)(full_model_pose);

		// show scores of start_pose and full_model_pose
		if ( option[ show_scores ]() ) {
			std::cout << "\n Score before building in missing residues:" << std::endl;
			scorefxn->show( std::cout, start_pose );
			std::cout << "\n Score after building in missing residues:" << std::endl;
			scorefxn->show( std::cout, full_model_pose );
		}

		// dump pdb of start_pose and full_model_pose
		if ( option[ dump_pdb ]() ) {
			std::string start_tag = tag + "_start";
			std::string full_model_tag = option[ virtualize_built ]() ? tag + "_full_model_virtual" : tag + "_full_model";
			start_pose.dump_pdb( start_tag + ".pdb" );
			full_model_pose.dump_pdb( full_model_tag + ".pdb" );
		}

		// write out to silent file
		stepwise::monte_carlo::output_to_silent_file( tag, silent_file_out, full_model_pose );

	}

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	build_full_model_test();
	protocols::viewer::clear_conformation_viewers();
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
		NEW_OPT( virtualize_built, "virtualize all newly built residues", false );
		NEW_OPT( show_scores, "show scores of pose before and after building in missing residues", false );
		NEW_OPT( dump_pdb, "dump pdb before and after building in missing residues", false );

		// setup
		core::init::init(argc, argv);
		option[ OptionKeys::chemical::patch_selectors ].push_back( "VIRTUAL_BASE" );
		option[ OptionKeys::chemical::patch_selectors ].push_back( "VIRTUAL_RIBOSE" );
		option[ OptionKeys::chemical::patch_selectors ].push_back( "VIRTUAL_RNA_RESIDUE" );
		option[ OptionKeys::chemical::patch_selectors ].push_back( "TERMINAL_PHOSPHATE" );

		// run
		protocols::viewer::viewer_main( my_main );
		exit( 0 );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
