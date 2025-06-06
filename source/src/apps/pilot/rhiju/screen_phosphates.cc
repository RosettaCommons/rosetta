// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file swa_monte_carlo.cc
/// @author Rhiju Das (rhiju@stanford.edu)

// libRosetta headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/FullModelPoseBuilder.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <protocols/rna/denovo/util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/rna/phosphate/MultiPhosphateSampler.hh>
#include <protocols/viewer/viewers.hh>

//////////////////////////////////////////////////
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh> // for option[ out::file::silent  ] and etc.
#include <basic/options/keys/in.OptionKeys.gen.hh> // for option[ in::file::tags ] and etc.
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>

//////////////////////////////////////////////////////////

// C++ headers
//#include <cstdlib>
#include <iostream>
#include <string>

using namespace core;
using namespace protocols;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using utility::vector1;

static basic::Tracer TR( "apps.pilot.rhiju.screen_phosphates" );

OPT_KEY( Boolean, screen_all_phosphates )
OPT_KEY( Boolean, do_screening )
OPT_KEY( IntegerVector, five_prime_phosphate_screen_res )
OPT_KEY( IntegerVector, three_prime_phosphate_screen_res )


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
screen_phosphates()
{
	using namespace core::pose;
	using namespace core::import_pose;
	using namespace core::scoring;
	using namespace core::id;
	using namespace core::chemical;
	using namespace core::pose::full_model_info;
	using namespace protocols::stepwise::modeler;
	using namespace protocols::stepwise::modeler::rna::phosphate;
	using namespace utility::file;

	// Following could be generalized to fa_standard, after recent unification, but
	// probably should wait for on-the-fly residue type generation.
	ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	PoseOP native_pose;
	if ( option[ in::file::native ].user() ) native_pose = get_pdb_and_cleanup( option[ in::file::native ](), rsd_set );

	// Following could go to a FullModelSetup class.
	// read starting pose(s) from disk
	utility::vector1< std::string > const & input_files = option[ in::file::s ]();
	utility::vector1< pose::PoseOP > input_poses;
	if ( input_files.size() == 0 ) input_poses.push_back( utility::pointer::make_shared< Pose >() ); // just a blank pose for now.
	for ( Size n = 1; n <= input_files.size(); n++ )  input_poses.push_back( get_pdb_and_cleanup( input_files[ n ], rsd_set ) );
	if ( option[ full_model::other_poses ].user() ) get_other_poses( input_poses, option[ full_model::other_poses ](), rsd_set );

	FullModelPoseBuilder builder;
	builder.set_input_poses( input_poses );
	builder.set_options( option );
	builder.initialize_further_from_options();
	builder.build(); // should update input_poses
	//fill_full_model_info_from_command_line( input_poses );  //FullModelInfo (minimal object needed for add/delete)

	// scorefunction
	core::scoring::ScoreFunctionOP scorefxn;
	if ( option[ score::weights ].user() ) scorefxn = get_score_function();
	else  scorefxn = ScoreFunctionFactory::create_score_function( "stepwise/rna/rna_res_level_energy.wts" );

	// actual pose to be sampled...
	PoseOP pose_op = input_poses[ 1 ];
	pose::Pose & pose = *pose_op;
	protocols::viewer::add_conformation_viewer ( pose.conformation(), "current", 500, 500 );
	protocols::stepwise::modeler::reroot_based_on_full_model_info( pose );

	fix_up_residue_type_variants( pose );

	// scorefxn->show( pose );

	Pose pose_copy = pose;

	MultiPhosphateSampler phosphate_sampler( pose );
	FullModelInfo const & full_model_info = const_full_model_info( pose );
	phosphate_sampler.set_scorefxn( scorefxn );
	phosphate_sampler.set_five_prime_phosphate_res(  full_model_info.full_to_sub( option[ five_prime_phosphate_screen_res  ]() ) );
	phosphate_sampler.set_three_prime_phosphate_res( full_model_info.full_to_sub( option[ three_prime_phosphate_screen_res ]() ) );
	phosphate_sampler.set_screen_all( option[ screen_all_phosphates ] );
	phosphate_sampler.sample_phosphates( pose_op );

	TR << "WITHOUT SCREEN PHOS (CONTROL) " << ( *scorefxn )( pose_copy ) << std::endl;
	// scorefxn->show( pose_copy );

	TR << "WITH SCREEN PHOS " << ( *scorefxn )( pose ) << std::endl;
	// scorefxn->show( pose );
	// protocols::rna::denovo::print_hbonds( pose );
	pose.dump_pdb( "SCREEN_POSE.pdb" );

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	clock_t const my_main_time_start( clock() );
	screen_phosphates();

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

		utility::vector1< core::Size > blank_size_vector;
		NEW_OPT( screen_all_phosphates, "Try screening phosphates at all residues where they could go", false );
		NEW_OPT( do_screening, "Actually do sample phosphates", false );
		NEW_OPT( five_prime_phosphate_screen_res, "Which residues (in full model numbering) to add 5' screenable phosphates to", blank_size_vector );
		NEW_OPT( three_prime_phosphate_screen_res, "Which residues (in full model numbering) to add 3' screenable phosphates to", blank_size_vector );

		option.add_relevant( in::file::fasta );
		option.add_relevant( in::file::input_res );
		option.add_relevant( in::file::native );
		option.add_relevant( out::file::silent );
		option.add_relevant( out::nstruct );
		option.add_relevant( score::weights );
		option.add_relevant( OptionKeys::full_model::extra_min_res );
		option.add_relevant( OptionKeys::full_model::sample_res );
		option.add_relevant( basic::options::OptionKeys::stepwise::monte_carlo::cycles );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::verbose_scores );
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
		option.add_relevant( OptionKeys::stepwise::rna::erraser );
		option.add_relevant( OptionKeys::full_model::rna::force_syn_chi_res_list );
		option.add_relevant( OptionKeys::stepwise::rna::virtual_sugar_keep_base_fixed );
		option.add_relevant( OptionKeys::stepwise::rna::force_centroid_interaction );
		option.add_relevant( OptionKeys::stepwise::num_random_samples );
		option.add_relevant( basic::options::OptionKeys::full_model::rna::bulge_res );
		option.add_relevant( basic::options::OptionKeys::full_model::rna::terminal_res );
		option.add_relevant( OptionKeys::rna::corrected_geo );

		devel::init(argc, argv);

		protocols::viewer::viewer_main( my_main );

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
}


