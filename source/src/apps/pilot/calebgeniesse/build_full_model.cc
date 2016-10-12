// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/calebgeniesse/build_full_model.cc
/// @brief Build in missing residues with stepwise::monte_carlo::build_full_model
/// @author Caleb Geniesse, geniesse@stanford.edu


// protocols
#include <protocols/farna/movers/RNA_LoopCloser.hh>
#include <protocols/farna/movers/RNA_LoopCloser.fwd.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/stepwise/monte_carlo/util.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/file_util.hh>
#include <protocols/stepwise/modeler/align/util.hh>

#include <protocols/farna/RNA_DeNovoProtocol.hh>
#include <protocols/farna/options/RNA_DeNovoProtocolOptions.hh>
#include <protocols/farna/util.hh>
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

using namespace core;
using namespace protocols;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::farna::options;
using namespace protocols::farna::setup;


OPT_KEY( Boolean, virtualize_built )
OPT_KEY( Boolean, show_scores )
OPT_KEY( Boolean, dump_pdb )
OPT_KEY( Boolean, refine_after_building )
OPT_KEY( Real, cluster_radius )

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.calebgeniesse.build_full_model" );



///////////////////////////////////////////////////////////////////////////////
bool
check_for_closeness( pose::Pose & pose_test, pose::Pose const & full_model_pose, Real const & rms_cutoff ) {
	Real const rms = stepwise::modeler::align::superimpose_with_stepwise_aligner( pose_test, full_model_pose, true /*superimpose_over_all_*/ );
	TR << TR.Cyan << "Calculated RMS to model: " << tag_from_pose( pose_test ) << " to be: " << rms << " and compared to " << rms_cutoff << TR.Reset << std::endl;
	return ( rms < rms_cutoff );
}

void
output_pose_list( utility::vector1< pose::PoseOP > pose_list, std::string const & silent_file ) {
	for ( Size i = 1; i <= pose_list.size(); ++i ) {
		pose::Pose const & p = *(pose_list[i]);
		std::string tag = tag_from_pose( p );
		stepwise::monte_carlo::output_to_silent_file( tag, silent_file, p);
	}
}



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
	using namespace protocols::farna;

	// setup residue types
	ResidueTypeSetCOP rsd_set;
	rsd_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	// setup score function
	core::scoring::ScoreFunctionOP scorefxn;
	if ( option[ score::weights ].user() ) {
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
	std::string const silent_file_rna_denovo_out = option[ out::file::silent ]() + ".rna_denovo.out";
	if ( option[ out::overwrite ]() ) {
		stepwise::modeler::remove_silent_file_if_it_exists( silent_file_rna_denovo_out );
	}
	std::string const silent_file_cluster_out = option[ out::file::silent ]() + ".cluster.out";
	if ( option[ out::overwrite ]() ) {
		stepwise::modeler::remove_silent_file_if_it_exists( silent_file_cluster_out );
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
	utility::vector1< PoseOP > pose_clusters;

	// clustering stuff..
	//Real cluster_rmsd = 0.0;
	//if ( option[ cluster_radius ].user() ) {
	// cluster_rmsd = option[ cluster_radius ];
	//}


	// setup movers
	farna::movers::RNA_LoopCloserOP rna_loop_closer( new farna::movers::RNA_LoopCloser() );
	rna_loop_closer->set_verbose( true );

	// setup RNA DeNovo stuff ...
	utility::vector1< core::pose::PoseOP > refine_pose_list;
	utility::vector1< std::string > chunk_silent_files;

	RNA_DeNovoProtocolOptionsOP rna_denovo_opts( new RNA_DeNovoProtocolOptions );
	rna_denovo_opts->initialize_from_command_line();
	rna_denovo_opts->set_chunk_silent_files( chunk_silent_files );
	rna_denovo_opts->set_silent_file( silent_file_rna_denovo_out );
	rna_denovo_opts->set_refine_pose( true );

	protocols::farna::RNA_DeNovoProtocol rna_de_novo_protocol( rna_denovo_opts );


	// iterate over input stream
	while ( input->has_another_pose() ) {

		// fill/score start_pose
		input->fill_pose( start_pose, *rsd_set );
		other_ops = nonconst_full_model_info( start_pose ).other_pose_list();
		stepwise::setup::fill_full_model_info_from_command_line( start_pose, other_ops );
		std::string const tag = tag_from_pose( start_pose );
		utility::vector1< Size > const full_model_cutpoint_open = const_full_model_info( start_pose ).cutpoint_open_in_full_model();

		for ( Size idx = 1; idx <= other_ops.size(); ++idx ) {
			if ( start_pose.total_residue() < other_ops[idx]->total_residue() ) {
				TR << "Switching focus to other pose at idx: " << TR.Cyan << idx << TR.Reset << "  (pose.nres=" << TR.Magenta << start_pose.total_residue() << TR.Reset << ", other.nres=" << TR.Magenta << other_ops[idx]->total_residue() << TR.Reset << ")" << std::endl;
				TR.Debug << "Start Pose N Residues: " << TR.Cyan << start_pose.total_residue() << TR.Reset << std::endl;
				TR.Debug << "Other Pose N Residues: " << TR.Cyan << other_ops[idx]->total_residue() << TR.Reset << std::endl;
				stepwise::modeler::switch_focus_to_other_pose( start_pose, idx /*scorefxn = 0 */);
				tag_into_pose( start_pose, tag );
				other_ops = nonconst_full_model_info( start_pose ).other_pose_list();
				break;
			}
		}

		nonconst_full_model_info( start_pose ).clear_other_pose_list();
		other_ops = nonconst_full_model_info( start_pose ).other_pose_list();

		FullModelInfo const & start_info = const_full_model_info( start_pose );
		(*scorefxn)(start_pose);

		// build full_model_pose
		if ( start_pose.total_residue() == start_info.full_sequence().size() ) {
			TR << "No missing residues in pose:  " << TR.Cyan << tag << TR.Reset << std::endl;
		} else {
			TR << "Building full model for pose: " << TR.Yellow << tag << TR.Reset << std::endl;
		}
		if ( other_ops.size() ) {
			TR << TR.Red << "[WARNING] other_pose_list.size() != 0, full_model_pose may be garbage!!! " << TR.Reset << std::endl;
		}

		// BUILD IN MISSING RESIDUES
		stepwise::monte_carlo::build_full_model( start_pose, full_model_pose );

		// refine via RNA DeNovo
		if ( option[ refine_after_building ]() ) {
			TR << "Refining full model for pose: " << TR.Cyan << tag << TR.Reset << std::endl;
			refine_pose_list.clear();
			refine_pose_list.push_back( PoseOP( new Pose( full_model_pose ) ) );
			rna_de_novo_protocol.set_refine_pose_list( refine_pose_list );
			rna_de_novo_protocol.apply( full_model_pose );
			stepwise::modeler::remove_silent_file_if_it_exists( silent_file_rna_denovo_out );
		}

		// compare full_model_info of start pose and full_model_pose
		FullModelInfo const & full_model_info = const_full_model_info( full_model_pose );
		//utility::vector1< Size > const & full_model_cutpoint_open = full_model_info.cutpoint_open_in_full_model();
		TR << "Cutpoint Open in full model:  " << TR.Magenta << full_model_cutpoint_open << TR.Reset << std::endl;

		TR << "Target Sequence: " << TR.Blue << start_info.full_sequence() << TR.Reset << std::endl;
		TR << "Start Sequence:  " << TR.Red << start_pose.sequence() << TR.Reset << std::endl;
		TR << "Built Sequence:  " << TR.Green << full_model_pose.sequence() << TR.Reset << std::endl;

		// virtualize built residues
		for ( Size i = 1; i <= full_model_info.res_list().size(); ++i ) {

			Size const & full_model_res = full_model_info.res_list()[i];
			Size other_idx = start_info.get_idx_for_other_pose_with_residue( full_model_res );

			if ( start_info.res_list().has_value( full_model_res ) || other_idx ) {

				// re-apply variants of residue in start_pose, to residue in full_model_pose
				ResidueOP start_rsd = 0;
				if ( other_idx ) {
					FullModelInfo const & other_info = const_full_model_info( *other_ops[other_idx] );
					start_rsd = ResidueOP( new Residue( other_ops[other_idx]->residue( other_info.full_to_sub( full_model_res ) ) ) );
				} else {
					start_rsd = ResidueOP( new Residue( start_pose.residue( start_info.full_to_sub( full_model_res ) ) ) );
				}
				utility::vector1< std::string > variant_types = start_rsd->type().properties().get_list_of_variants();

				for ( Size j = 1; j <= variant_types.size(); ++j ) {
					VariantType const & variant_type = start_rsd->type().properties().get_variant_from_string( variant_types[j] );
					add_variant_type_to_pose_residue( full_model_pose, variant_type, full_model_res );
				}

			} else {

				// close cutpoints
				other_idx = start_info.get_idx_for_other_pose_with_residue( full_model_res + 1 );
				/*if ( start_info.res_list().has_value( full_model_res + 1 ) || other_idx ) {
				TR << "Closing loop at residue " << TR.Yellow << full_model_res << TR.Reset << std::endl;
				add_variant_type_to_pose_residue( full_model_pose, CUTPOINT_LOWER, full_model_res );
				add_variant_type_to_pose_residue( full_model_pose, CUTPOINT_UPPER, full_model_res + 1 );
				rna_loop_closer->apply( full_model_pose ); //, full_model_res );
				}*/

				// apply virtual_rna_residue variant to residue in full_model_pose
				if ( option[ virtualize_built ]() ) {
					remove_variant_type_from_pose_residue( full_model_pose, BLOCK_STACK_BELOW, full_model_res );
					remove_variant_type_from_pose_residue( full_model_pose, BLOCK_STACK_ABOVE, full_model_res );
					remove_variant_type_from_pose_residue( full_model_pose, CUTPOINT_LOWER, full_model_res );
					if ( full_model_res < full_model_pose.total_residue() ) {
						remove_variant_type_from_pose_residue( full_model_pose, CUTPOINT_UPPER, full_model_res + 1 );
					}
					if ( full_model_res > 1 && !full_model_cutpoint_open.has_value( full_model_res ) ) {
						TR << "Virtualizing residue " << TR.Yellow << full_model_res << TR.Reset << std::endl;
						apply_virtual_rna_residue_variant_type( full_model_pose, full_model_res, false /* apply_check */ );
					} else {
						TR << "Cutpoint Open at residue " << TR.Magenta << full_model_res << TR.Reset << std::endl;
					}

				}
			}
		}

		// update pdb info/ score full_model_pose
		update_pdb_info_from_full_model_info( full_model_pose );

		// show scores of start_pose and full_model_pose
		if ( option[ show_scores ]() ) {
			(*scorefxn)(full_model_pose);
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


	if ( pose_clusters.size() > 0 ) {
		output_pose_list( pose_clusters, silent_file_cluster_out );
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
		NEW_OPT( refine_after_building, "refine pose via RNA DeNovo after building in missing residues", false );
		NEW_OPT( cluster_radius, "cluster radius (clustering off by default)", 0.0 );

		// setup
		core::init::init(argc, argv);
		option[ OptionKeys::chemical::patch_selectors ].push_back( "VIRTUAL_BASE" );
		option[ OptionKeys::chemical::patch_selectors ].push_back( "VIRTUAL_RIBOSE" );
		option[ OptionKeys::chemical::patch_selectors ].push_back( "VIRTUAL_RNA_RESIDUE" );
		option[ OptionKeys::chemical::patch_selectors ].push_back( "TERMINAL_PHOSPHATE" );
		option[ OptionKeys::chemical::patch_selectors ].push_back( "LOWER_TERMINUS_VARIANT" );

		// run
		protocols::viewer::viewer_main( my_main );
		exit( 0 );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
