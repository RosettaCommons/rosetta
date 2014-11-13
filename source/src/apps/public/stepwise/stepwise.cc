// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file swa_monte_carlo.cc
/// @author Rhiju Das (rhiju@stanford.edu)

// libRosetta headers
#include <core/types.hh>
#include <core/chemical/util.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/init/init.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/util.hh>
#include <core/io/rna/RNA_DataReader.cc> // temporary, for scoring RNA chemical mapping data. Move into core?
#include <protocols/stepwise/full_model_info/FullModelInfoSetupFromCommandLine.hh>
#include <protocols/stepwise/monte_carlo/StepWiseMonteCarlo.hh>
#include <protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.hh>
#include <protocols/stepwise/monte_carlo/util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/file_util.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/viewer/viewers.hh>

//////////////////////////////////////////////////
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>
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

using namespace protocols;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using utility::vector1;

static thread_local basic::Tracer TR( "apps.pilot.rhiju.stepwise_monte_carlo" );


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
stepwise_monte_carlo()
{
  using namespace core::pose;
  using namespace core::scoring;
  using namespace core::chemical;
  using namespace core::pose::full_model_info;
  using namespace protocols::stepwise;
  using namespace protocols::stepwise::modeler;
	using namespace protocols::stepwise::full_model_info;
  using namespace protocols::stepwise::monte_carlo;
  using namespace utility::file;

	// Following could be generalized to fa_standard, after recent unification, but
	// probably should wait for on-the-fly residue type generation.
	ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	// scorefunction
	core::scoring::ScoreFunctionOP scorefxn;
	if ( option[ score::weights ].user() ) scorefxn = get_score_function();
	else  scorefxn = ScoreFunctionFactory::create_score_function( "stepwise/rna/rna_res_level_energy.wts" );
	if ( option[ OptionKeys::constraints::cst_file ].user() && !scorefxn->has_nonzero_weight( atom_pair_constraint ) ) scorefxn->set_weight( atom_pair_constraint, 1.0 );

	PoseOP native_pose, align_pose;
	initialize_native_and_align_pose( native_pose, align_pose, rsd_set );
	PoseOP pose_op = initialize_pose_and_other_poses_from_command_line( rsd_set );
	pose::Pose & pose = *pose_op;
	// temporary, for scoring RNA chemical mapping data. Move into initalize_pose?
	core::io::rna::get_rna_data_info( pose, option[ basic::options::OptionKeys::rna::data_file ](), scorefxn );

	// Get rid of this commented code when it is incorporated into a unit test.
	//	test_merge_and_slice_with_two_helix_test_case( input_poses, scorefxn ); exit( 0 );

	// actual pose to be sampled...
	if ( pose.total_residue() > 0 ) ( *scorefxn )( pose );
	Vector center_vector = ( align_pose != 0 ) ? get_center_of_mass( *align_pose ) : Vector( 0.0 );
	protocols::viewer::add_conformation_viewer ( pose.conformation(), "current", 500, 500, false, ( align_pose != 0 ), center_vector );

	StepWiseMonteCarlo stepwise_monte_carlo( scorefxn );
	protocols::stepwise::monte_carlo::options::StepWiseMonteCarloOptionsOP options( new protocols::stepwise::monte_carlo::options::StepWiseMonteCarloOptions );
	bool const do_preminimize_move = option[ OptionKeys::stepwise::preminimize ]();
	bool const test_move = option[ OptionKeys::stepwise::move ].user() || do_preminimize_move;
	if ( test_move ) options->set_output_minimized_pose_list( true );
	options->initialize_from_command_line();
	stepwise_monte_carlo.set_options( options );
	stepwise_monte_carlo.set_native_pose( align_pose ); //allows for alignment to be to non-native
	stepwise_monte_carlo.set_move( SWA_Move( option[ OptionKeys::stepwise::move ](), const_full_model_info( pose ).full_model_parameters() ) );
	stepwise_monte_carlo.set_enumerate( option[ OptionKeys::stepwise::enumerate ]());
	stepwise_monte_carlo.set_do_preminimize_move( do_preminimize_move );

	std::string const silent_file = option[ out::file::silent ]();
	if ( option[ out::overwrite ]() ) remove_silent_file_if_it_exists( silent_file );
	stepwise_monte_carlo.set_out_path( FileName( silent_file ).path() );

	std::string out_tag;
	Pose start_pose = pose;

	// main loop
	for ( Size n = 1; n <= Size( option[ out::nstruct ]() ); n++ ) {
		if ( !get_out_tag( out_tag, n, silent_file ) ) continue;
		TR << std::endl << TR.Green << "Embarking on structure " << n << " of " << option[ out::nstruct ]() << TR.Reset << std::endl;
		pose = start_pose;
		stepwise_monte_carlo.set_model_tag( out_tag );
 		stepwise_monte_carlo.apply( pose );
		if (!options->output_minimized_pose_list()) output_to_silent_file( out_tag, silent_file, pose, native_pose, option[ OptionKeys::stepwise::superimpose_over_all ](), true /*rms_fill*/ );
	}

	if ( do_preminimize_move ) pose.dump_pdb( "PREPACK.pdb" );
}

///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	clock_t const my_main_time_start( clock() );
	stepwise_monte_carlo();
	protocols::viewer::clear_conformation_viewers();
	std::cout << "Total time to run " << static_cast<Real>( clock() - my_main_time_start ) / CLOCKS_PER_SEC << " seconds." << std::endl;
  exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
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
		option.add_relevant( OptionKeys::full_model::extra_min_res );
		option.add_relevant( OptionKeys::full_model::jump_res );
		option.add_relevant( OptionKeys::full_model::root_res );
		option.add_relevant( OptionKeys::full_model::virtual_sugar_res );
		option.add_relevant( OptionKeys::full_model::cutpoint_open );
		option.add_relevant( OptionKeys::full_model::cutpoint_closed );
		option.add_relevant( OptionKeys::full_model::sample_res );
		option.add_relevant( OptionKeys::stepwise::superimpose_over_all );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::cycles );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::skip_deletions );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::add_delete_frequency );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::minimize_single_res_frequency );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::switch_focus_frequency );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::just_min_after_mutation_frequency );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::allow_internal_hinge_moves );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::allow_internal_local_moves );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::allow_skip_bulge );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::temperature );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::allow_variable_bond_geometry );
		option.add_relevant( OptionKeys::stepwise::move );
		option.add_relevant( OptionKeys::stepwise::num_random_samples );
		option.add_relevant( OptionKeys::stepwise::num_pose_minimize );
		option.add_relevant( OptionKeys::stepwise::align_pdb );
		option.add_relevant( OptionKeys::stepwise::enumerate );
		option.add_relevant( OptionKeys::stepwise::preminimize );
		option.add_relevant( OptionKeys::stepwise::rna::erraser );
		option.add_relevant( OptionKeys::stepwise::rna::force_centroid_interaction );
		option.add_relevant( OptionKeys::stepwise::rna::rebuild_bulge_mode );
		option.add_relevant( OptionKeys::full_model::rna::force_syn_chi_res_list );
		option.add_relevant( OptionKeys::full_model::rna::bulge_res );
		option.add_relevant( OptionKeys::full_model::rna::terminal_res );
		option.add_relevant( OptionKeys::stepwise::atr_rep_screen );
		option.add_relevant( OptionKeys::stepwise::protein::allow_virtual_side_chains );
		option.add_relevant( OptionKeys::rna::corrected_geo );
		option.add_relevant( OptionKeys::rna::data_file );

		core::init::init(argc, argv);

		option[ OptionKeys::chemical::patch_selectors ].push_back( "VIRTUAL_SIDE_CHAIN" ); // for protein side-chain packing/virtualization
		option[ OptionKeys::chemical::patch_selectors ].push_back( "VIRTUAL_RIBOSE" ); // for skip-nucleotide moves.
		option[ OptionKeys::chemical::patch_selectors ].push_back( "VIRTUAL_BASE" ); // for chemical mapping & bulge.
		option[ OptionKeys::chemical::patch_selectors ].push_back( "VIRTUAL_RNA_RESIDUE" ); // for bulge
		option[ OptionKeys::chemical::patch_selectors ].push_back( "PEPTIDE_CAP" ); // N_acetylated.txt and C_methylamidated.txt
		option[ OptionKeys::chemical::patch_selectors ].push_back( "TERMINAL_PHOSPHATE" ); // 5prime_phosphate and 3prime_phosphate

		protocols::viewer::viewer_main( my_main );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}



