// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>

#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/OptionKeyList.hh>

#include <basic/Tracer.hh>

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::stepwise::modeler;

static basic::Tracer TR( "protocols.stepwise.monte_carlo.options.StepWiseMonteCarloOptions" );

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace options {

//Constructor
StepWiseMonteCarloOptions::StepWiseMonteCarloOptions():
	verbose_scores_( false ),
	integration_test_mode_( false ),
	force_centroid_interaction_( true ),
	sampler_max_centroid_distance_( 0.0 ),
	use_phenix_geo_( true ),
	skip_deletions_( false ),
	erraser_( true ),
	cycles_( 500 ),
	add_proposal_density_factor_( 1.0 ),
	minimize_single_res_frequency_( 0.0 ),
	just_min_after_mutation_frequency_( 0.5 ),
	temperature_( 1.0 ),
	allow_split_off_( true ),
	virtual_sugar_keep_base_fixed_( false ),
	virtual_sugar_do_minimize_( false /*true*/ ),
	make_movie_( false ),
	sampler_perform_phosphate_pack_( true ),
	force_phosphate_instantiation_( true ),
	virtualize_packable_moieties_in_screening_pose_( true ),
	rebuild_bulge_mode_( false ),
	tether_jump_( true ),
	local_redock_only_( true ),
	skip_coord_constraints_( false ),
	filter_native_big_bins_( false ),
	allow_virtual_o2prime_hydrogens_( false ),
	allow_virtual_side_chains_( false ),
	n_sample_( 18 ),
	protein_prepack_( false ),
	o2prime_legacy_mode_( false ),
	recover_low_( false ),
	enumerate_( false ),
	preminimize_( false ),
	skip_preminimize_( false ),
	new_move_selector_( false ),
	test_all_moves_( false ),
	save_times_( false ),
	use_precomputed_library_( true ),
	minimize_after_delete_( true ),
	use_first_jump_for_submotif_( false ),
	designing_with_noncanonicals_( false ),
	checkpoint_( false ),
	checkpointing_frequency_( 0 ),
	continue_until_none_missing_( false ),
	eval_base_pairs_( false ),
	superimpose_over_all_( false ),
	force_moving_res_for_erraser_( false )
{
	StepWiseBasicOptions::initialize_variables();
	set_silent_file( "default.out" );
}

//Destructor
StepWiseMonteCarloOptions::~StepWiseMonteCarloOptions()
{}

/// @brief copy constructor
StepWiseMonteCarloOptions::StepWiseMonteCarloOptions( StepWiseMonteCarloOptions const & src ) :
	ResourceOptions( src ),
	StepWiseBasicOptions( src ),
	StepWiseMoveSelectorOptions( src )
{
	*this = src;
}

/// @brief clone the options
StepWiseMonteCarloOptionsOP
StepWiseMonteCarloOptions::clone() const
{
	return StepWiseMonteCarloOptionsOP( new StepWiseMonteCarloOptions( *this ) );
}

void
StepWiseMonteCarloOptions::initialize_from_options_collection( utility::options::OptionCollection const & options ) {

	StepWiseBasicOptions::initialize_from_options_collection( options );
	StepWiseMoveSelectorOptions::initialize_from_options_collection( options );

	set_silent_file( options[ OptionKeys::out::file::silent ]() );
	set_verbose_scores( options[ OptionKeys::stepwise::monte_carlo::verbose_scores ]() );
	set_integration_test_mode( options[ OptionKeys::stepwise::rna::integration_test ]() );
	set_skip_deletions( options[ OptionKeys::stepwise::monte_carlo::skip_deletions ]() );
	set_num_pose_minimize( options[ OptionKeys::stepwise::num_pose_minimize ]() );
	set_erraser( options[ OptionKeys::stepwise::rna::erraser ]() );
	set_cycles( options[ OptionKeys::stepwise::monte_carlo::cycles ]() );
	set_add_proposal_density_factor( options[ OptionKeys::stepwise::monte_carlo::add_proposal_density_factor ]() );
	set_minimize_single_res_frequency( options[ OptionKeys::stepwise::monte_carlo::minimize_single_res_frequency ]() );
	set_just_min_after_mutation_frequency( options[ OptionKeys::stepwise::monte_carlo::just_min_after_mutation_frequency ]() );

	set_allow_split_off( options[ OptionKeys::stepwise::monte_carlo::allow_split_off ]() );
	set_temperature( options[ OptionKeys::stepwise::monte_carlo::temperature ]() );
	set_virtual_sugar_keep_base_fixed( options[ OptionKeys::stepwise::rna::virtual_sugar_keep_base_fixed ]() );
	if ( options[ OptionKeys::stepwise::rna::virtual_sugar_do_minimize ].user() ) {
		set_virtual_sugar_do_minimize( options[ OptionKeys::stepwise::rna::virtual_sugar_do_minimize ]() );
	}
	force_centroid_interaction_ = true; // note default is different from stepwise enumeration
	if ( options[ OptionKeys::stepwise::rna::force_centroid_interaction ].user() ) set_force_centroid_interaction( options[ OptionKeys::stepwise::rna::force_centroid_interaction ]() );
	sampler_max_centroid_distance_ = options[ OptionKeys::stepwise::rna::sampler_max_centroid_distance ]();
	set_use_phenix_geo(  options[ OptionKeys::rna::corrected_geo ] );
	set_make_movie( options[ OptionKeys::stepwise::monte_carlo::make_movie ]() );
	sampler_perform_phosphate_pack_ = options[ OptionKeys::stepwise::rna::sampler_perform_phosphate_pack ]();
	force_phosphate_instantiation_ = options[ OptionKeys::stepwise::rna::force_phosphate_instantiation ]();
	rebuild_bulge_mode_ = options[ OptionKeys::stepwise::rna::rebuild_bulge_mode ]();
	tether_jump_ = options[ OptionKeys::stepwise::rna::tether_jump ]();
	o2prime_legacy_mode_ = options[ OptionKeys::stepwise::rna::o2prime_legacy_mode ]();
	recover_low_ = options[ OptionKeys::stepwise::monte_carlo::recover_low ]();
	enumerate_ = options[ OptionKeys::stepwise::enumerate ]();
	preminimize_ = options[ OptionKeys::stepwise::preminimize ]();
	skip_preminimize_ = options[ OptionKeys::stepwise::skip_preminimize ]();
	new_move_selector_ = options[ OptionKeys::stepwise::new_move_selector ]();
	test_all_moves_ = options[ OptionKeys::stepwise::test_all_moves ]();
	save_times_ = options[ OptionKeys::out::save_times ]();
	use_precomputed_library_ = options[ OptionKeys::stepwise::monte_carlo::use_precomputed_library ]();
	local_redock_only_ = options[ OptionKeys::stepwise::monte_carlo::local_redock_only ]();
	skip_coord_constraints_ = options[ OptionKeys::stepwise::protein::skip_coord_constraints ]();
	filter_native_big_bins_ = options[ OptionKeys::stepwise::protein::filter_native_big_bins ]();
	allow_virtual_o2prime_hydrogens_ = options[ OptionKeys::stepwise::rna::allow_virtual_o2prime_hydrogens ]();
	allow_virtual_side_chains_ = options[ OptionKeys::stepwise::protein::allow_virtual_side_chains ]();
	n_sample_ = options[ OptionKeys::stepwise::protein::n_sample ]();
	protein_prepack_ = options[ OptionKeys::stepwise::protein::protein_prepack ]();
	virtualize_packable_moieties_in_screening_pose_ = options[ OptionKeys::stepwise::virtualize_packable_moieties_in_screening_pose ]();
	use_first_jump_for_submotif_ = options[ OptionKeys::stepwise::monte_carlo::use_first_jump_for_submotif ]();
	exclude_submotifs_ = options[ OptionKeys::stepwise::monte_carlo::exclude_submotifs ]();
	designing_with_noncanonicals_ = options[ OptionKeys::stepwise::monte_carlo::designing_with_noncanonicals ]();
	checkpoint_ = options[ OptionKeys::stepwise::monte_carlo::checkpointing_frequency ]() != 0;
	checkpointing_frequency_ = options[ OptionKeys::stepwise::monte_carlo::checkpointing_frequency ]();
	continue_until_none_missing_ = options[ OptionKeys::stepwise::monte_carlo::continue_until_none_missing ]();
	eval_base_pairs_ = options[ basic::options::OptionKeys::rna::evaluate_base_pairs ]();
	superimpose_over_all_ = options[ basic::options::OptionKeys::stepwise::superimpose_over_all ]();
	force_moving_res_for_erraser_ = options[ basic::options::OptionKeys::stepwise::force_moving_res_for_erraser ]();

	if ( test_all_moves_ ) {
		set_num_random_samples( 0 );
		set_minimize_after_delete( false );
	}
}

void
StepWiseMonteCarloOptions::list_options_read( utility::options::OptionKeyList & opts ) {

	StepWiseBasicOptions::list_options_read( opts );
	StepWiseMoveSelectorOptions::list_options_read( opts );

	opts + OptionKeys::out::file::silent
		+ OptionKeys::stepwise::monte_carlo::verbose_scores
		+ OptionKeys::stepwise::rna::integration_test
		+ OptionKeys::stepwise::monte_carlo::skip_deletions
		+ OptionKeys::stepwise::num_pose_minimize
		+ OptionKeys::stepwise::rna::erraser
		+ OptionKeys::stepwise::monte_carlo::cycles
		+ OptionKeys::stepwise::monte_carlo::add_proposal_density_factor
		+ OptionKeys::stepwise::monte_carlo::minimize_single_res_frequency
		+ OptionKeys::stepwise::monte_carlo::just_min_after_mutation_frequency
		+ OptionKeys::stepwise::monte_carlo::allow_split_off
		+ OptionKeys::stepwise::monte_carlo::temperature
		+ OptionKeys::stepwise::rna::virtual_sugar_keep_base_fixed
		+ OptionKeys::stepwise::rna::virtual_sugar_do_minimize
		+ OptionKeys::stepwise::rna::force_centroid_interaction
		+ OptionKeys::stepwise::rna::sampler_max_centroid_distance
		+ OptionKeys::rna::corrected_geo
		+ OptionKeys::stepwise::monte_carlo::make_movie
		+ OptionKeys::stepwise::rna::sampler_perform_phosphate_pack
		+ OptionKeys::stepwise::rna::force_phosphate_instantiation
		+ OptionKeys::stepwise::rna::rebuild_bulge_mode
		+ OptionKeys::stepwise::rna::tether_jump
		+ OptionKeys::stepwise::rna::o2prime_legacy_mode
		+ OptionKeys::stepwise::monte_carlo::recover_low
		+ OptionKeys::stepwise::enumerate
		+ OptionKeys::stepwise::preminimize
		+ OptionKeys::stepwise::skip_preminimize
		+ OptionKeys::stepwise::new_move_selector
		+ OptionKeys::stepwise::test_all_moves
		+ OptionKeys::out::save_times
		+ OptionKeys::stepwise::monte_carlo::use_precomputed_library
		+ OptionKeys::stepwise::monte_carlo::local_redock_only
		+ OptionKeys::stepwise::protein::skip_coord_constraints
		+ OptionKeys::stepwise::protein::filter_native_big_bins
		+ OptionKeys::stepwise::rna::allow_virtual_o2prime_hydrogens
		+ OptionKeys::stepwise::protein::allow_virtual_side_chains
		+ OptionKeys::stepwise::protein::n_sample
		+ OptionKeys::stepwise::protein::protein_prepack
		+ OptionKeys::stepwise::virtualize_packable_moieties_in_screening_pose
		+ OptionKeys::stepwise::monte_carlo::use_first_jump_for_submotif
		+ OptionKeys::stepwise::monte_carlo::exclude_submotifs
		+ OptionKeys::stepwise::monte_carlo::designing_with_noncanonicals
		+ OptionKeys::stepwise::monte_carlo::checkpointing_frequency
		+ OptionKeys::stepwise::monte_carlo::continue_until_none_missing
		+ OptionKeys::rna::evaluate_base_pairs
		+ OptionKeys::stepwise::superimpose_over_all
		+ OptionKeys::stepwise::force_moving_res_for_erraser;
}

///////////////////////////////////////////////////////////////////
void
StepWiseMonteCarloOptions::initialize_from_command_line() {
	initialize_from_options_collection( option );
}


//////////////////////////////////////////////////////////////////////////////////////
// Wait, would be smarter to just take advantage of clone(), right?
protocols::stepwise::modeler::options::StepWiseModelerOptionsOP
StepWiseMonteCarloOptions::setup_modeler_options() const{

	using namespace protocols::stepwise::options;
	using namespace protocols::stepwise::modeler::options;
	StepWiseModelerOptionsOP options( new StepWiseModelerOptions );

	// copy over StepWiseBasicOptions -- both StepWiseMonteCarloOptions and StepWiseModelerOptions inherit from it.
	StepWiseBasicOptions & options_basic_new( *options );
	StepWiseBasicOptions const & options_basic_this( *this );
	options_basic_new = options_basic_this;

	// general
	options->set_num_pose_minimize( 0 );
	options->set_num_random_samples( num_random_samples() );
	if ( num_pose_minimize() > 0 ) options->set_num_pose_minimize( num_pose_minimize() );
	options->set_rmsd_screen( rmsd_screen() );
	options->set_VDW_rep_screen_info( VDW_rep_screen_info() ); // caleb when you see a conflict, just get rid of this line.
	options->set_atr_rep_screen( atr_rep_screen() );

	// Might be smarter to inherit the following options from the same options parent class,
	// instead of setting manually here, where we could forget something.
	options->set_choose_random( !enumerate_ );
	options->set_virtualize_packable_moieties_in_screening_pose( virtualize_packable_moieties_in_screening_pose() );
	if ( preminimize_ ) options->set_use_packer_instead_of_rotamer_trials( true );  // to match SWA for proteins.
	if ( enumerate_ && !erraser() ) {
		options->set_output_minimized_pose_list( true ); // to match legacy SWA (turn off for erraser)
	}
	// protein-specific
	options->set_skip_coord_constraints( skip_coord_constraints() );
	options->set_filter_native_big_bins( filter_native_big_bins() );
	options->set_allow_virtual_o2prime_hydrogens( allow_virtual_o2prime_hydrogens() ); // RNA
	options->set_allow_virtual_side_chains( allow_virtual_side_chains() ); // protein
	options->set_n_sample( n_sample() );
	options->set_prepack( protein_prepack() );
	options->set_expand_loop_takeoff( true );

	// rna-specific
	options->set_integration_test_mode( integration_test_mode() );
	options->set_force_centroid_interaction( force_centroid_interaction() );
	options->set_sampler_max_centroid_distance( sampler_max_centroid_distance() );
	options->set_use_phenix_geo( use_phenix_geo() );
	options->set_kic_modeler_if_relevant( erraser() );
	options->set_virtual_sugar_keep_base_fixed( virtual_sugar_keep_base_fixed() );
	options->set_virtual_sugar_do_minimize( virtual_sugar_do_minimize() );
	if ( rebuild_bulge_mode_ ) options->set_force_centroid_interaction( false );
	options->set_tether_jump( tether_jump() );
	options->set_o2prime_legacy_mode( o2prime_legacy_mode() );
	options->set_sampler_perform_phosphate_pack( sampler_perform_phosphate_pack() );
	options->set_force_phosphate_instantiation( force_phosphate_instantiation() );
	options->set_allow_rebuild_bulge_mode( false );
	options->set_virtual_sugar_do_screens( false );

	return options;
}


} //options
} //monte_carlo
} //stepwise
} //protocols
