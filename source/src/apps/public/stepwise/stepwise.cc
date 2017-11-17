// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file stepwise.cc
/// @author Rhiju Das (rhiju@stanford.edu)
/// @author Andrew Watkins (amw579@stanford.edu)

// libRosetta headers
#include <core/types.hh>
#include <core/chemical/util.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/rna/RNA_DataReader.hh> // temporary, for scoring RNA chemical mapping data. Move into core?
#include <core/io/silent/util.hh>
#include <protocols/scoring/VDW_CachedRepScreenInfo.hh>
#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>
#include <protocols/rna/setup/RNA_JobDistributor.hh>
#include <protocols/rna/setup/RNA_CSA_JobDistributor.hh>
#include <protocols/rna/setup/RNA_MonteCarloJobDistributor.hh>
#include <protocols/stepwise/monte_carlo/StepWiseMonteCarlo.hh>
#include <protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.hh>
#include <protocols/stepwise/monte_carlo/util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/monte_carlo/submotif/SubMotifLibrary.hh>
#include <protocols/viewer/viewers.hh>

#include <protocols/jd3/full_model/FullModelJobQueen.hh>
#include <protocols/jd3/full_model/MoverAndFullModelJob.hh>
#include <protocols/jd3/full_model_inputters/PDBFullModelInputter.hh>
#include <protocols/jd3/JobDistributor.hh>
#include <protocols/jd3/JobDistributorFactory.hh>
#include <protocols/jd3/Job.fwd.hh>
#include <protocols/jd3/LarvalJob.fwd.hh>
#include <protocols/jd3/JobResult.fwd.hh>

//////////////////////////////////////////////////
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh> // for option[ in::file::tags ] and etc.
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh> // for option[ out::file::silent  ] and etc.
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>

#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>

//////////////////////////////////////////////////////////
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <list>

using namespace core;
using namespace pose;
using namespace protocols;
using namespace protocols::jd3;
using namespace protocols::jd3::full_model;
using namespace protocols::jd3::full_model_inputters;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace basic::options::OptionKeys::stepwise::monte_carlo;
using utility::vector1;

static basic::Tracer TR( "apps.public.stepwise.stepwise" );

OPT_KEY( Boolean, use_legacy_stepwise_job_distributor )

core::scoring::ScoreFunctionOP
get_stepwise_score_function( OptionCollection const & option ) {
	using namespace core::scoring;
	ScoreFunctionOP scorefxn;
	if ( option[ score::weights ].user() ) {
		scorefxn = get_score_function();
	} else if ( option[ OptionKeys::stepwise::lores ]() && !option[ OptionKeys::rna::denovo::minimize_rna ]() ) {
		scorefxn = ScoreFunctionFactory::create_score_function( "stepwise/rna/rna_lores_for_stepwise.wts" );
	} else {
		// force user to make a conscious choice.
		utility_exit_with_message( "Please specify scorefunction.\n\n  Current release for RNA     : -score:weights stepwise/rna/rna_res_level_energy7beta.wts\n  Last published for RNA      : -score:weights stepwise/rna/rna_res_level_energy4.wts -restore_talaris_behavior\n  Current test for RNA/protein: -score:weights stepwise/stepwise_res_level_energy.wts\n" );
	}
	if ( option[ OptionKeys::constraints::cst_file ].user() && !scorefxn->has_nonzero_weight( atom_pair_constraint ) ) {
		scorefxn->set_weight( atom_pair_constraint, 1.0 );
	}
	return scorefxn;
}

class StepWiseJobQueen : public FullModelJobQueen {

public:
	StepWiseJobQueen():
		FullModelJobQueen()
	{
		utility::options::OptionKeyList opts;
		//protocols::stepwise::options::StepWiseBasicOptions::list_options_read( opts );
		//protocols::stepwise::options::StepWiseMoveSelectorOptions::list_options_read( opts );
		protocols::stepwise::monte_carlo::options::StepWiseMonteCarloOptions::list_options_read( opts );
		PDBFullModelInputter::list_options_read( opts );
		core::scoring::list_read_options_in_get_score_function( opts );
		add_options( opts );

		//add_option( basic::options::OptionKeys::score::weights );
		add_option( basic::options::OptionKeys::stepwise::monte_carlo::csa::csa_bank_size );
		add_option( basic::options::OptionKeys::stepwise::lores );
		add_option( basic::options::OptionKeys::rna::denovo::minimize_rna );
		add_option( basic::options::OptionKeys::constraints::cst_file );
		add_option( basic::options::OptionKeys::stepwise::move );
		add_option( basic::options::OptionKeys::out::file::silent );
		add_option( basic::options::OptionKeys::out::nstruct );
		add_option( basic::options::OptionKeys::out::overwrite );
		add_option( basic::options::OptionKeys::stepwise::superimpose_over_all );
	}

	virtual
	JobOP
	complete_larval_job_maturation(
		protocols::jd3::LarvalJobCOP larval_job,
		utility::options::OptionCollectionCOP job_options,
		utility::vector1< JobResultCOP > const & // input_job_results
	) {
		using namespace core::scoring;
		using namespace core::pose;
		using namespace core::chemical;
		using namespace core::pose::full_model_info;
		using namespace core::import_pose;
		using namespace protocols::stepwise;
		using namespace protocols::stepwise::setup;
		using namespace protocols::stepwise::modeler;
		using namespace protocols::stepwise::monte_carlo::submotif;
		using namespace protocols::stepwise::monte_carlo;
		using namespace protocols::stepwise::monte_carlo::mover;
		using namespace protocols::stepwise::monte_carlo::options;

		MoverAndFullModelJobOP mature_job( new MoverAndFullModelJob );

		core::pose::PoseOP pose = pose_for_job( larval_job, *job_options );
		mature_job->pose( pose );

		utility::options::OptionCollection const & jo = *job_options;

		ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

		// AMW: pass this scorefunction to the JQ
		ScoreFunctionOP scorefxn = get_stepwise_score_function( jo );
		if ( jo[ score::weights ].user() ) scorefxn = get_score_function();

		// If we care about the 'viewability' of the pose, center it -- obviously this
		// isn't super-useful unless we have a BIG starting pose where residue 1 is
		// likely far from the COM.
#ifdef GL_GRAPHICS
		pose->center();
#endif
		PoseOP native_pose, align_pose;
		initialize_native_and_align_pose( native_pose, align_pose, rsd_set, pose );
		// temporary, for scoring RNA chemical mapping data. Move into initalize_pose?
		core::io::rna::get_rna_data_info( *pose, option[ basic::options::OptionKeys::rna::data_file ](), scorefxn );

		// setup test move specified via -stepwise:move option
		StepWiseMove const test_move( jo[ OptionKeys::stepwise::move ](), const_full_model_info( *pose ).full_model_parameters() );

		// actual pose to be sampled... do not score pose is user has specified a move
		if ( pose->size() > 0 && test_move.move_type() == NO_MOVE ) ( *scorefxn )( *pose );
		Vector center_vector = ( align_pose != nullptr ) ? get_center_of_mass( *align_pose ) : Vector( 0.0 );
		protocols::viewer::add_conformation_viewer ( pose->conformation(), "current", 500, 500, false, ( align_pose != 0 ), center_vector );

		StepWiseMonteCarloOP stepwise_monte_carlo( new StepWiseMonteCarlo( scorefxn ) );
		stepwise_monte_carlo->set_align_pose( align_pose ); //allows for alignment to be to non-native
		stepwise_monte_carlo->set_move( test_move );

		StepWiseMonteCarloOptionsOP options( new StepWiseMonteCarloOptions );
		options->initialize_from_options_collection( *job_options );
		stepwise_monte_carlo->set_options( options );
		if ( ( options->from_scratch_frequency() > 0.0 || const_full_model_info( *pose ).other_pose_list().size() > 0 ) && !scorefxn->has_nonzero_weight( other_pose ) ) scorefxn->set_weight( other_pose, 1.0 ); // critical if more than one pose shows up and focus switches...

		std::string const silent_file = jo[ out::file::silent ]();
		if ( jo[ out::overwrite ]() ) core::io::silent::remove_silent_file_if_it_exists( silent_file );
		stepwise_monte_carlo->set_out_path( utility::file::FileName( silent_file ).path() );
		stepwise_monte_carlo->set_submotif_library( SubMotifLibraryCOP( new SubMotifLibrary( rsd_set, options->lores() /*include_submotifs_from_jump_library*/, options->use_first_jump_for_submotif(), options->exclude_submotifs() ) ) );

		// used to just be the stepwise job distributor that got this.
		// did we have behavior that relied on these two having a different
		// sense? yes... if align was null we made it -s! THIS native was never
		// meant for rms fill calculations before.
		stepwise_monte_carlo->set_native_pose( native_pose );

		mature_job->pose( pose );

		// It takes a mover and pose, but the pose must have a FullModelInfo.
		// Don't worry about the latter check for now.
		mature_job->mover( stepwise_monte_carlo );

		return mature_job;
	}
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
stepwise_monte_carlo()
{
	// AMW: THIS will determine whether it's MPI vs serial etc.
	// Elsewhere, in the DAG setup, we will decide MC vs CSA
	protocols::jd3::JobDistributorOP jd = protocols::jd3::JobDistributorFactory::create_job_distributor();
	protocols::jd3::JobQueenOP queen( new StepWiseJobQueen );
	jd->go( queen );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
stepwise_monte_carlo_legacy()
{
	using namespace core::pose;
	using namespace core::scoring;
	using namespace core::chemical;
	using namespace core::pose::full_model_info;
	using namespace core::import_pose;
	using namespace protocols::stepwise;
	using namespace protocols::stepwise::modeler;
	using namespace protocols::stepwise::monte_carlo::submotif;
	using namespace protocols::stepwise::setup;
	using namespace protocols::rna::setup;
	using namespace protocols::stepwise::monte_carlo;
	using namespace protocols::stepwise::monte_carlo::mover;
	using namespace protocols::stepwise::monte_carlo::options;
	using namespace utility::file;

	ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	ScoreFunctionOP scorefxn = get_stepwise_score_function( option );

	PoseOP native_pose, align_pose;
	PoseOP pose_op = initialize_pose_and_other_poses_from_command_line( rsd_set );
	protocols::scoring::fill_vdw_cached_rep_screen_info_from_command_line( *pose_op );
	pose::Pose & pose = *pose_op;

	// If we care about the 'viewability' of the pose, center it -- obviously this
	// isn't super-useful unless we have a BIG starting pose where residue 1 is
	// likely far from the COM.
#ifdef GL_GRAPHICS
	pose.center();
#endif

	initialize_native_and_align_pose( native_pose, align_pose, rsd_set, pose_op );
	// temporary, for scoring RNA chemical mapping data. Move into initalize_pose?
	core::io::rna::get_rna_data_info( pose, option[ basic::options::OptionKeys::rna::data_file ](), scorefxn );

	// Get rid of this commented code when it is incorporated into a unit test.
	// test_merge_and_slice_with_two_helix_test_case( input_poses, scorefxn ); exit( 0 );

	// setup test move specified via -stepwise:move option
	StepWiseMove const test_move( option[ OptionKeys::stepwise::move ](), const_full_model_info( pose ).full_model_parameters() );

	// actual pose to be sampled... do not score pose is user has specified a move
	if ( pose.size() > 0 && test_move.move_type() == NO_MOVE ) ( *scorefxn )( pose );
	Vector center_vector = ( align_pose != 0 ) ? get_center_of_mass( *align_pose ) : Vector( 0.0 );
	protocols::viewer::add_conformation_viewer ( pose.conformation(), "current", 500, 500, false, ( align_pose != 0 ), center_vector );

	StepWiseMonteCarloOP stepwise_monte_carlo( new StepWiseMonteCarlo( scorefxn ) );
	stepwise_monte_carlo->set_native_pose( align_pose ); //allows for alignment to be to non-native
	stepwise_monte_carlo->set_move( test_move );

	StepWiseMonteCarloOptionsOP options( new StepWiseMonteCarloOptions );
	options->initialize_from_command_line();
	stepwise_monte_carlo->set_options( options );
	if ( ( options->from_scratch_frequency() > 0.0 || const_full_model_info( *pose_op ).other_pose_list().size() > 0 ) && !scorefxn->has_nonzero_weight( other_pose ) ) {
		scorefxn->set_weight( other_pose, 1.0 ); // critical if more than one pose shows up and focus switches...
	}

	std::string const silent_file = option[ out::file::silent ]();
	if ( option[ out::overwrite ]() ) {
		core::io::silent::remove_silent_file_if_it_exists( silent_file );
	}
	stepwise_monte_carlo->set_out_path( FileName( silent_file ).path() );
	stepwise_monte_carlo->set_submotif_library( SubMotifLibraryCOP( new SubMotifLibrary( rsd_set, options->lores() /*include_submotifs_from_jump_library*/, options->use_first_jump_for_submotif(), options->exclude_submotifs() ) ) );

	// main loop
	RNA_JobDistributorOP stepwise_job_distributor( new RNA_MonteCarloJobDistributor( stepwise_monte_carlo, silent_file, option[ out::nstruct ]() ) );
	if ( option[ csa::csa_bank_size ].user() ) {
		stepwise_job_distributor = RNA_JobDistributorOP( new RNA_CSA_JobDistributor( stepwise_monte_carlo, silent_file, option[ out::nstruct ](), option[ csa::csa_bank_size ](), option[ csa::csa_rmsd ](), option[ csa::csa_output_rounds ](), option[ csa::annealing ]() ) );
	}
	stepwise_job_distributor->set_native_pose( native_pose );
	stepwise_job_distributor->set_superimpose_over_all( option[ OptionKeys::stepwise::superimpose_over_all ]() );
	stepwise_job_distributor->initialize( pose );

	while ( stepwise_job_distributor->has_another_job() ) {
		stepwise_job_distributor->apply( pose );
	}
}

///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	clock_t const my_main_time_start( clock() );
	if ( option[ use_legacy_stepwise_job_distributor ].value() ) {
		stepwise_monte_carlo_legacy();
	} else {
		// Check that the user is not supplying an -out:file:silent with a path
		std::string out_file_silent = option[ out::file::silent ].value();
		if ( out_file_silent.find( '/' ) != std::string::npos ) {
			utility_exit_with_message( "Error: you are trying to specify -out:file:silent in a way that specifies a path, which is not compatible with the new job distributor. Either use -out:file:path as well as -out:file:silent, or pass -use_legacy_stepwise_job_distributor." );
		}
		stepwise_monte_carlo();
	}
	protocols::viewer::clear_conformation_viewers();
	std::cout << "Total time to run " << static_cast<Real>( clock() - my_main_time_start ) / CLOCKS_PER_SEC << " seconds." << std::endl;
	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		NEW_OPT( use_legacy_stepwise_job_distributor, "Use the legacy RNA_MonteCarloJobDistributor", true );

		std::cout << std::endl << "Basic usage:  " << argv[0] << "  -fasta <fasta file with sequence> -s <start pdb> -input_res <input pdb1> [ -native <native pdb file> ] " << std::endl;
		std::cout << std::endl << " Type -help for full slate of options." << std::endl << std::endl;

		option.add_relevant( in::file::fasta );
		option.add_relevant( in::file::input_res );
		option.add_relevant( in::file::native );
		option.add_relevant( out::file::silent );
		option.add_relevant( out::overwrite );
		option.add_relevant( out::nstruct );
		option.add_relevant( score::weights );
		option.add_relevant( constraints::cst_file );
		option.add_relevant( OptionKeys::full_model::other_poses );
		option.add_relevant( OptionKeys::full_model::jump_res );
		option.add_relevant( OptionKeys::full_model::root_res );
		option.add_relevant( OptionKeys::full_model::virtual_sugar_res );
		option.add_relevant( OptionKeys::full_model::cutpoint_open );
		option.add_relevant( OptionKeys::full_model::cutpoint_closed );
		option.add_relevant( OptionKeys::full_model::sample_res );
		option.add_relevant( OptionKeys::full_model::motif_mode );
		option.add_relevant( OptionKeys::full_model::rna::force_syn_chi_res_list );
		option.add_relevant( OptionKeys::full_model::rna::force_anti_chi_res_list );
		option.add_relevant( OptionKeys::full_model::rna::force_south_sugar_list );
		option.add_relevant( OptionKeys::full_model::rna::force_north_sugar_list );
		option.add_relevant( OptionKeys::full_model::rna::bulge_res );
		option.add_relevant( OptionKeys::full_model::extra_min_res );
		option.add_relevant( OptionKeys::full_model::rna::terminal_res );
		option.add_relevant( OptionKeys::full_model::rna::block_stack_above_res );
		option.add_relevant( OptionKeys::full_model::rna::block_stack_below_res );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::cycles );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::skip_deletions );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::add_delete_frequency );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::minimize_single_res_frequency );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::switch_focus_frequency );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::just_min_after_mutation_frequency );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::submotif_frequency );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::allow_internal_hinge_moves );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::allow_internal_local_moves );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::allow_skip_bulge );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::temperature );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::make_movie );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::allow_variable_bond_geometry );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::csa::csa_bank_size );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::csa::csa_rmsd );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::csa::csa_output_rounds );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::csa::annealing );
		option.add_relevant( OptionKeys::stepwise::superimpose_over_all );
		option.add_relevant( OptionKeys::stepwise::move );
		option.add_relevant( OptionKeys::stepwise::num_random_samples );
		option.add_relevant( OptionKeys::stepwise::num_pose_minimize );
		option.add_relevant( OptionKeys::stepwise::align_pdb );
		option.add_relevant( OptionKeys::stepwise::enumerate );
		option.add_relevant( OptionKeys::stepwise::preminimize );
		option.add_relevant( OptionKeys::stepwise::atr_rep_screen );
		option.add_relevant( OptionKeys::stepwise::min_type );
		option.add_relevant( OptionKeys::stepwise::min_tolerance );
		option.add_relevant( OptionKeys::stepwise::virtualize_free_moieties_in_native );
		option.add_relevant( OptionKeys::stepwise::new_move_selector );
		option.add_relevant( OptionKeys::stepwise::rna::erraser );
		option.add_relevant( OptionKeys::stepwise::rna::force_centroid_interaction );
		option.add_relevant( OptionKeys::stepwise::rna::rebuild_bulge_mode );
		option.add_relevant( OptionKeys::stepwise::rna::integration_test );
		option.add_relevant( OptionKeys::stepwise::protein::allow_virtual_side_chains );
		option.add_relevant( OptionKeys::rna::corrected_geo );
		option.add_relevant( OptionKeys::rna::data_file );

		devel::init(argc, argv);

		if ( option[ OptionKeys::stepwise::superimpose_over_all ].user() ) {
			std::cout << "The use of -superimpose_over_all is deprecated. The behavior in question now defaults to TRUE and is turned off by providing a particular residue that is part of an anchoring input domain as -alignment_anchor_res." << std::endl;
		}

		protocols::viewer::viewer_main( my_main );

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}


