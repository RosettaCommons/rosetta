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
#include <protocols/viewer/viewers.hh>
#include <protocols/stepwise/monte_carlo/util.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/file_util.hh>
#include <protocols/stepwise/modeler/align/util.hh>

#include <protocols/farna/FARNA_Optimizer.hh>
#include <protocols/farna/RNA_DeNovoProtocol.hh>
#include <protocols/farna/options/RNA_DeNovoProtocolOptions.hh>
#include <protocols/farna/util.hh>
#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>


// core
#include <core/types.hh>
#include <devel/init.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/sequence/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
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
#include <core/scoring/Energies.hh>
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
using namespace core::scoring;
using namespace core::chemical;
using namespace core::conformation;
using namespace core::io::silent;
using namespace core::import_pose::pose_stream;
using namespace core::pose;
using namespace core::pose::rna;
using namespace core::pose::full_model_info;
using namespace protocols::farna;

OPT_KEY( Boolean, virtualize_built )
OPT_KEY( Boolean, show_scores )
OPT_KEY( Boolean, dump_pdb )
OPT_KEY( Boolean, refine_after_building )
OPT_KEY( Boolean, fragment_assembly_mode )
OPT_KEY( Boolean, caleb_legacy )
OPT_KEY( Boolean, focus_on_largest_pose )

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.calebgeniesse.build_full_model" );

void
output_seq_with_residue_emphasis( std::string const & seq, Size const res ) {
	TR.Trace << "Seq: ";
	for ( Size ii = 1; ii <= seq.size(); ++ii ) {
		if ( ii == res ) {
			TR.Trace << TR.Blue << seq[ ii - 1 ] << TR.Reset;
		} else {
			TR.Trace << seq[ ii - 1 ];
		}
	}
	TR.Trace << std::endl;
}

class BuildFullModel {
public:

	BuildFullModel();

	void
	make_built_residues_repulsive(
		Pose const & start_pose,
		Pose & full_model_pose,
		utility::vector1< PoseOP > const & //other_ops
	);

	void
	make_built_residues_virtual(
		Pose const & start_pose,
		Pose & full_model_pose,
		utility::vector1< PoseOP > const & other_ops
	);

	void
	fill_and_sample_full_model(
		Pose & start_pose, // because of show()?
		utility::vector1< PoseOP > const & other_ops
	);

private:
	protocols::farna::RNA_DeNovoProtocolOP rna_de_novo_protocol_ = nullptr;
	core::pose::PoseCOP native_pose_ = nullptr;
};

BuildFullModel::BuildFullModel() {

	// 1. Setup RNA_DeNovoProtocol

	std::string const silent_file_rna_denovo_out = option[ out::file::silent ]() + ".rna_denovo.out";
	if ( option[ out::overwrite ]() ) {
		stepwise::modeler::remove_silent_file_if_it_exists( silent_file_rna_denovo_out );
	}
	utility::vector1< std::string > chunk_silent_files;

	RNA_DeNovoProtocolOptionsOP rna_denovo_opts( new RNA_DeNovoProtocolOptions );
	rna_denovo_opts->initialize_from_command_line();
	rna_denovo_opts->set_chunk_silent_files( chunk_silent_files );
	rna_denovo_opts->set_silent_file( silent_file_rna_denovo_out );
	rna_denovo_opts->set_refine_pose( option[ caleb_legacy ]() );
	if ( !option[ caleb_legacy ]() ) {
		// Temp - debugging. Set to ~5000 production.
		rna_denovo_opts->set_monte_carlo_cycles( 500 );
		rna_denovo_opts->set_user_defined_cycles( true );
	}
	rna_de_novo_protocol_ = protocols::farna::RNA_DeNovoProtocolOP( new protocols::farna::RNA_DeNovoProtocol( rna_denovo_opts ) );

	// 2. Setup native pose
	if ( option[ in::file::native ].user() ) {
		ResidueTypeSetCOP rsd_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		native_pose_ = stepwise::setup::get_pdb_with_full_model_info( option[ in::file::native ], rsd_set ); //import_pose::pose_from_file( option[ in::file::native ] );
	}
}

void
BuildFullModel::make_built_residues_repulsive(
	Pose const & start_pose,
	Pose & full_model_pose,
	utility::vector1< PoseOP > const & //other_ops
) {
	FullModelInfo const & start_info = const_full_model_info( start_pose );
	FullModelInfo const & full_model_info = const_full_model_info( full_model_pose );
	utility::vector1< Size > const & full_model_cutpoint_open = const_full_model_info( start_pose ).cutpoint_open_in_full_model();

	std::string const & full_sequence = start_info.full_sequence();

	// Add rep-only patch to built residues.
	// virtualize built residues
	for ( Size const full_model_res : full_model_info.res_list() ) {

		Size other_idx = start_info.get_idx_for_other_pose_with_residue( full_model_res );

		if ( start_info.res_list().has_value( full_model_res ) || other_idx ) {
			/*
			// TODO: check domain map for input stuff (nonzero?).

			// re-apply variants of residue in start_pose, to residue in full_model_pose
			ResidueOP start_rsd = 0;
			if ( other_idx ) {
			FullModelInfo const & other_info = const_full_model_info( *other_ops[other_idx] );
			start_rsd = ResidueOP( new Residue( other_ops[other_idx]->residue( other_info.full_to_sub( full_model_res ) ) ) );
			} else {
			start_rsd = ResidueOP( new Residue( start_pose.residue( start_info.full_to_sub( full_model_res ) ) ) );
			}
			utility::vector1< std::string > variant_types = start_rsd->type().properties().get_list_of_variants();

			// It's not safe to give residues their old variants back. What if they imply
			// a connectivity state that isn't true anymore?

			//for ( Size j = 1; j <= variant_types.size(); ++j ) {
			// VariantType const & variant_type = start_rsd->type().properties().get_variant_from_string( variant_types[j] );
			// add_variant_type_to_pose_residue( full_model_pose, variant_type, full_model_info.full_to_sub( full_model_res ) );
			//}
			*/
		} else {
			if ( full_model_res > 1 && !full_model_cutpoint_open.has_value( full_model_res ) ) {
				output_seq_with_residue_emphasis( full_sequence, full_model_res );
				TR.Trace << "Repulsive-only-ing residue " << TR.Yellow << full_model_res << TR.Reset << std::endl;

				// Add the REPLONLY variant using the same method as
				// add_virtual_rna_residue_type
				pose::add_variant_type_to_pose_residue( full_model_pose, REPLONLY, full_model_res );

				// Check in terms of transitions between input domains
				// Ensure not a cutpoint open in full_model_info

				bool is_cutpoint_closed = false;
				if ( full_model_pose.residue_type( full_model_res ).has_variant_type( CUTPOINT_LOWER ) ) {
					runtime_assert( full_model_pose.residue_type( full_model_res + 1 ).has_variant_type( CUTPOINT_UPPER ) );
					is_cutpoint_closed = true;
				}

				bool const is_cutpoint_open = ( full_model_pose.fold_tree().is_cutpoint( full_model_res ) && !is_cutpoint_closed );

				if ( ( full_model_pose.size() != full_model_res ) &&  ( !is_cutpoint_open ) ) {
					pose::add_variant_type_to_pose_residue( full_model_pose, REPL_PHOSPHATE, full_model_res + 1 );
				}

			} else {
				TR << "Cutpoint Open at residue " << TR.Magenta << full_model_res << TR.Reset << std::endl;
			}
		}
	}
}

void BuildFullModel::make_built_residues_virtual(
	Pose const & start_pose,
	Pose & full_model_pose,
	utility::vector1< PoseOP > const & other_ops
) {
	FullModelInfo const & start_info = const_full_model_info( start_pose );
	FullModelInfo const & full_model_info = const_full_model_info( full_model_pose );
	utility::vector1< Size > const & full_model_cutpoint_open = const_full_model_info( start_pose ).cutpoint_open_in_full_model();

	for ( Size const full_model_res : full_model_info.res_list() ) {

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
			//other_idx = start_info.get_idx_for_other_pose_with_residue( full_model_res + 1 );

			// apply virtual_rna_residue variant to residue in full_model_pose
			if ( option[ virtualize_built ]() ) {
				remove_variant_type_from_pose_residue( full_model_pose, BLOCK_STACK_BELOW, full_model_info.full_to_sub(full_model_res) );
				remove_variant_type_from_pose_residue( full_model_pose, BLOCK_STACK_ABOVE, full_model_info.full_to_sub(full_model_res) );
				remove_variant_type_from_pose_residue( full_model_pose, CUTPOINT_LOWER, full_model_info.full_to_sub(full_model_res) );
				if ( full_model_res < full_model_pose.total_residue() ) {
					remove_variant_type_from_pose_residue( full_model_pose, CUTPOINT_UPPER, full_model_info.full_to_sub(full_model_res + 1) );

					// Maybe now that we have removed the cutpoint we should add a bond.
					full_model_pose.conformation().set_polymeric_connection(full_model_info.full_to_sub(full_model_res), full_model_info.full_to_sub(full_model_res+1));
				}

				if ( full_model_res > 1 && !full_model_cutpoint_open.has_value( full_model_info.full_to_sub(full_model_res) ) ) {
					TR << "Virtualizing residue " << TR.Yellow << full_model_info.full_to_sub(full_model_res) << TR.Reset << std::endl;
					apply_virtual_rna_residue_variant_type( full_model_pose, full_model_info.full_to_sub(full_model_res), false /* apply_check */ );
				} else {
					TR << "Cutpoint Open at residue " << TR.Magenta << full_model_info.full_to_sub(full_model_res) << TR.Reset << std::endl;
				}
			}
		}
	}
}

void
BuildFullModel::fill_and_sample_full_model(
	Pose & start_pose, // because of show() and score
	utility::vector1< PoseOP > const & other_ops
) {
	// setup score function
	core::scoring::ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( "stepwise/rna/rna_res_level_energy4" );

	std::string const tag = tag_from_pose( start_pose );

	(*scorefxn)(start_pose);

	FullModelInfo const & start_info = const_full_model_info( start_pose );
	core::scoring::EnergiesOP start_energies = start_pose.energies().clone();

	bool will_need_fragment_assembly_stage = false;
	if ( start_pose.total_residue() == start_info.full_sequence().size() ) {
		TR << "No missing residues in pose:  " << TR.Cyan << tag << TR.Reset << std::endl;
	} else {
		TR << "Building full model for pose: " << TR.Yellow << tag << TR.Reset << std::endl;
		will_need_fragment_assembly_stage = true;
	}

	FullModelParametersOP starting_full_model_parameters = const_full_model_info( start_pose ).full_model_parameters()->clone();

	utility::vector1< Size > domain = core::pose::full_model_info::figure_out_pose_domain_map_const( start_pose );

	// Graphics
	//protocols::viewer::add_conformation_viewer( full_model_pose.conformation(), "current", 500, 500, false, false );
	// Temporarily off -- will set up e.g. 2000 windows

	// BUILD IN MISSING RESIDUES

	Pose full_model_pose;
	stepwise::monte_carlo::build_full_model( start_pose, full_model_pose );

	if ( option[ fragment_assembly_mode ]() && ! option[ caleb_legacy ]() ) {

		make_built_residues_repulsive( start_pose, full_model_pose, other_ops );

		// Use fragment assembly to refine the positions of
		if ( will_need_fragment_assembly_stage ) {
			// Be extra-sure that the INPUT_DOMAIN and FIXED_DOMAIN are up to date.
			FullModelParametersOP fmp = nonconst_full_model_info( full_model_pose ).full_model_parameters()->clone();

			fmp->set_parameter( INPUT_DOMAIN, domain );
			fmp->set_parameter( FIXED_DOMAIN, domain );
			nonconst_full_model_info( full_model_pose ).set_full_model_parameters( fmp );

			utility::vector1< PoseOP > pose_list;
			pose_list.emplace_back( new Pose( full_model_pose ) );

			// TODO: consider scoring function decisions.
			static ScoreFunctionOP sfxn( new ScoreFunction );
			sfxn->set_weight( rna_vdw, 1.0 );
			sfxn->set_weight( linear_chainbreak, 1.0 );

			protocols::farna::FARNA_Optimizer farna_opt( pose_list, sfxn, 5000 );
			farna_opt.apply( full_model_pose );
		}
	}

	// refine via RNA_DeNovoProtocol -- old method, also did minimization
	if ( option[ refine_after_building ]() && option[ caleb_legacy ] ) {
		TR << "Refining full model for pose: " << TR.Cyan << tag << TR.Reset << std::endl;
		utility::vector1< core::pose::PoseOP > refine_pose_list;
		refine_pose_list.emplace_back( new Pose( full_model_pose ) );
		rna_de_novo_protocol_->set_refine_pose_list( refine_pose_list );
		rna_de_novo_protocol_->apply( full_model_pose );
		stepwise::modeler::remove_silent_file_if_it_exists( rna_de_novo_protocol_->options()->silent_file() );
	}

	//utility::vector1< Size > const & full_model_cutpoint_open = full_model_info.cutpoint_open_in_full_model();
	utility::vector1< Size > const & full_model_cutpoint_open = const_full_model_info( start_pose ).cutpoint_open_in_full_model();
	if ( !full_model_cutpoint_open.empty() ) {
		TR << "Cutpoint Open in full model:  " << TR.Magenta << full_model_cutpoint_open << TR.Reset << std::endl;
	}
	TR << "Target Sequence: " << TR.Blue << start_info.full_sequence() << TR.Reset << std::endl;
	TR << "Start Sequence:  " << TR.Red << start_pose.sequence() << TR.Reset << std::endl;
	TR << "Built Sequence:  " << TR.Green << full_model_pose.sequence() << TR.Reset << std::endl;

	// virtualize built residues
	if ( ! option[ fragment_assembly_mode ]() && option[ caleb_legacy ]() ) {
		make_built_residues_virtual( start_pose, full_model_pose, other_ops );
	}

	// update pdb info/ score full_model_pose
	update_pdb_info_from_full_model_info( full_model_pose );

	// AMW: We explicitly don't want these scores. Rather, we want to re-insert
	// the old scores (to follow).
	//(*scorefxn)(full_model_pose);

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

	// We need to put back in the old FullModelParameters -- because we have
	// job-specific data in it now (the input domain etc.)
	// and that's not desired -- you'd like that to be a constant over the
	// whole run...
	nonconst_full_model_info( full_model_pose ).set_full_model_parameters( starting_full_model_parameters );

	// We need to explicitly pop in the starting poses's energies for a couple reasons:
	// 1. The new FullModelPose may have some energy issues (repulsive residues)
	// ...
	full_model_pose.set_new_energies_object( start_energies );

	// write out to silent file
	// Note that we re-call 'build_full_model' to get rms_fill in this function below.
	// This is fine, since we have no missing residues at this point. It's a
	// no-op.
	stepwise::monte_carlo::output_to_silent_file( tag, option[ out::file::silent ](), full_model_pose, native_pose_, false, true );

	protocols::viewer::clear_conformation_viewers();
}

void
build_full_model_test()
{
	// setup residue types
	ResidueTypeSetCOP rsd_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	// setup silent_files
	if ( !option[ in::file::silent ].user() ) {
		utility_exit_with_message( "Error: Must supply a silent file! Use -in::file::silent <silent file>" );
	}

	if ( option[ out::overwrite ]() ) {
		stepwise::modeler::remove_silent_file_if_it_exists( option[ out::file::silent ]() );
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

	BuildFullModel bfm;

	// iterate over input stream
	while ( input->has_another_pose() ) {

		Pose start_pose;

		input->fill_pose( start_pose, *rsd_set );

		utility::vector1< PoseOP > other_ops = nonconst_full_model_info( start_pose ).other_pose_list();
		stepwise::setup::fill_full_model_info_from_command_line( start_pose, other_ops );

		std::string const tag = tag_from_pose( start_pose );

		// This particular legacy behavior forced a switch in focus to the
		// largest pose. This was important because we deleted all other residues.
		if ( option[ focus_on_largest_pose ]() || option[ caleb_legacy ]() ) {
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
		}

		if ( option[ caleb_legacy ]() ) {
			nonconst_full_model_info( start_pose ).clear_other_pose_list();
			other_ops = nonconst_full_model_info( start_pose ).other_pose_list();
			if ( !other_ops.empty() ) {
				TR << TR.Red << "[WARNING] other_pose_list not empty, full_model_pose may be garbage!!! " << TR.Reset << std::endl;
			}
		}

		bfm.fill_and_sample_full_model( start_pose, other_ops );
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
		NEW_OPT( fragment_assembly_mode, "Fill in missing residues using fragment assembly", false );
		NEW_OPT( caleb_legacy, "Make Caleb's legacy decision to nuke OTHER_POSEs", false );
		NEW_OPT( focus_on_largest_pose, "Switch focus to the largest pose.", false );

		// setup
		devel::init(argc, argv);
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
