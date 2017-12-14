// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/denovo/options/RNA_FragmentMonteCarloOptions.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/rna/denovo/options/RNA_FragmentMonteCarloOptions.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh> // for superimpose_over_all
#include <basic/database/open.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/OptionKeyList.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.rna.denovo.options.RNA_FragmentMonteCarloOptions" );

using namespace basic::options;
using namespace basic::options::OptionKeys;

namespace protocols {
namespace rna {
namespace denovo {
namespace options {

//Constructor
RNA_FragmentMonteCarloOptions::RNA_FragmentMonteCarloOptions():
	rounds_( 10 ),
	monte_carlo_cycles_( 0 ), /* will be reset later */
	user_defined_cycles_( false ), /* will change to true if set_monte_carlo_cycles() is called */
	minimize_structure_( false ),
	relax_structure_( false ),
	temperature_( 2.0 ),
	ignore_secstruct_( false ),
	chainbreak_weight_( -1.0 ), /* use rna/denovo/rna_lores.wts number unless user specified. -1.0 is never really used. */
	linear_chainbreak_weight_( -1.0 ), /* use rna/denovo/rna_lores.wts number unless user specified. -1.0 is never really used. */
	close_loops_( true ),
	close_loops_after_each_move_( false ),
	allow_bulge_( false ),
	allow_consecutive_bulges_( false ),
	jump_change_frequency_( 0.1 ),
	autofilter_( false ),
	autofilter_score_quantile_( 0.20 ),
	titrate_stack_bonus_( true ),
	root_at_first_rigid_body_( false ),
	dock_each_chunk_( false ),
	dock_each_chunk_per_chain_( false ),
	center_jumps_in_single_stranded_( false ),
	suppress_bp_constraint_( 1.0 ),
	filter_lores_base_pairs_( false ),
	filter_lores_base_pairs_early_( false ),
	filter_chain_closure_( true ),
	filter_chain_closure_distance_( 6.0 ), /* in Angstroms. This is pretty loose!*/
	filter_chain_closure_halfway_( true ),
	staged_constraints_( false ),
	output_score_frequency_( 0 ),
	output_jump_res_( 0 ),
	output_jump_o3p_to_o5p_( false ),
	output_jump_chainbreak_( false ),
	output_rotation_vector_( false ),
	output_jump_target_xyz_( 0.0 ),
	output_jump_reference_RT_string_( "" ),
	save_jump_histogram_( false ),
	jump_histogram_boxsize_( 0.0 ),
	jump_histogram_binwidth_( 0.0 ),
	jump_histogram_binwidth_rotvector_( 0.0 ),
	filter_vdw_( false ),
	vdw_rep_screen_include_sidechains_( false ),
	gradual_constraints_( true ),
	ramp_rnp_vdw_( false ),
	grid_vdw_weight_( 1.0 ),
	convert_protein_centroid_( true ),
	rna_protein_docking_( false ),
	small_docking_moves_( false ),
	docking_move_size_( 1.0 ),
	rna_protein_docking_legacy_( false ),
	rna_protein_docking_freq_( 0.4 ),
	docking_( false ),
	randomize_init_rnp_( true ),
	rnp_high_res_relax_( true ),
	rnp_high_res_cycles_( 10 ),
	rnp_pack_first_( false ),
	rnp_ramp_rep_( false ),
	rnp_min_first_( false ),
	dock_into_density_legacy_( false ),
	new_fold_tree_initializer_( false ),
	initial_structures_provided_( false ),
	simple_rmsd_cutoff_relax_( false ),
	refine_pose_( false ),
	override_refine_pose_rounds_( false ),
	refine_native_get_good_FT_( false ),
	bps_moves_( false ),
	disallow_bps_at_extra_min_res_( false ),
	allow_fragment_moves_in_bps_( false ),
	frag_size_( 0 ),
	// following is odd, but note that core::scoring::rna::chemical_shift machinery also checks global options system.
	use_chem_shift_data_( option[ OptionKeys::score::rna::rna_chemical_shift_exp_data].user() ),
	superimpose_over_all_( false ),
	fixed_stems_( false ),
	all_rna_fragments_file_( basic::database::full_name("sampling/rna/RICHARDSON_RNA09.torsions") ),
	rna_params_file_( "" ),
	jump_library_file_( basic::database::full_name("sampling/rna/1jj2_RNA_jump_library.dat" ) ),
	rmsd_screen_( 0.0 ),
	disallow_realign_( true )
{}

//Destructor
RNA_FragmentMonteCarloOptions::~RNA_FragmentMonteCarloOptions() = default;

/// @brief copy constructor
RNA_FragmentMonteCarloOptions::RNA_FragmentMonteCarloOptions( RNA_FragmentMonteCarloOptions const & src ) :
	ResourceOptions( src ),
	RNA_BasicOptions( src ),
	RNA_MinimizerOptions( src )
{
	*this = src;
}

/// @brief clone the options
RNA_FragmentMonteCarloOptionsOP
RNA_FragmentMonteCarloOptions::clone() const
{
	return RNA_FragmentMonteCarloOptionsOP( new RNA_FragmentMonteCarloOptions( *this ) );
}

///////////////////////////////////////////////////////////////////
void
RNA_FragmentMonteCarloOptions::initialize_from_command_line() {
	initialize_from_options( basic::options::option );
}

void
RNA_FragmentMonteCarloOptions::initialize_from_options( utility::options::OptionCollection const & opts ) {
	RNA_MinimizerOptions::initialize_from_options( opts ); // includes RNA_BasicOptions

	if ( opts[ basic::options::OptionKeys::rna::denovo::cycles ].user() ) {
		set_monte_carlo_cycles( opts[ basic::options::OptionKeys::rna::denovo::cycles ]() );
		set_user_defined_cycles( true );
	}
	set_rounds( opts[ basic::options::OptionKeys::rna::denovo::rounds ]() );
	minimize_structure_ = opts[ basic::options::OptionKeys::rna::denovo::minimize_rna ]();
	relax_structure_ = opts[ basic::options::OptionKeys::rna::denovo::relax_rna ]();
	allow_bulge_ = opts[ basic::options::OptionKeys::rna::denovo::allow_bulge ]();
	set_temperature( opts[ basic::options::OptionKeys::rna::denovo::temperature ] );
	set_ignore_secstruct( opts[ basic::options::OptionKeys::rna::denovo::ignore_secstruct ] );
	set_jump_change_frequency( opts[ basic::options::OptionKeys::rna::denovo::jump_change_frequency ] );
	set_close_loops( opts[ basic::options::OptionKeys::rna::denovo::close_loops] );
	set_close_loops_after_each_move( opts[ basic::options::OptionKeys::rna::denovo::close_loops_after_each_move ] );
	set_staged_constraints( opts[ basic::options::OptionKeys::rna::denovo::staged_constraints ] ) ;
	set_filter_lores_base_pairs( opts[ basic::options::OptionKeys::rna::denovo::filter_lores_base_pairs] );
	set_filter_lores_base_pairs_early( opts[ basic::options::OptionKeys::rna::denovo::filter_lores_base_pairs_early] );
	set_suppress_bp_constraint( opts[ basic::options::OptionKeys::rna::denovo::suppress_bp_constraint ]() );
	set_filter_chain_closure( opts[ basic::options::OptionKeys::rna::denovo::filter_chain_closure ]() );
	set_filter_chain_closure_distance( opts[ basic::options::OptionKeys::rna::denovo::filter_chain_closure_distance ]() );
	set_filter_chain_closure_halfway( opts[ basic::options::OptionKeys::rna::denovo::filter_chain_closure_halfway ]() );
	set_vary_bond_geometry( opts[ basic::options::OptionKeys::rna::vary_geometry ] );
	set_root_at_first_rigid_body( opts[ basic::options::OptionKeys::rna::denovo::root_at_first_rigid_body ] );
	set_dock_each_chunk( opts[ basic::options::OptionKeys::rna::denovo::dock_each_chunk ]() );
	set_dock_each_chunk_per_chain( opts[ basic::options::OptionKeys::rna::denovo::dock_each_chunk_per_chain ]() );
	set_center_jumps_in_single_stranded( opts[ basic::options::OptionKeys::rna::denovo::center_jumps_in_single_stranded ]() );
	set_autofilter( opts[ basic::options::OptionKeys::rna::denovo::autofilter ] );
	set_bps_moves( opts[ basic::options::OptionKeys::rna::denovo::bps_moves ] );
	set_disallow_bps_at_extra_min_res( opts[ basic::options::OptionKeys::rna::denovo::disallow_bps_at_extra_min_res ] );
	set_allow_fragment_moves_in_bps( opts[ basic::options::OptionKeys::rna::denovo::allow_fragment_moves_in_bps ] );
	set_frag_size( opts[ basic::options::OptionKeys::rna::denovo::frag_size ]() );

	set_output_score_frequency( opts[ basic::options::OptionKeys::rna::denovo::out::output_score_frequency ]() );
	set_output_score_file( opts[ basic::options::OptionKeys::rna::denovo::out::output_score_file ]() );
	set_output_jump_res( opts[ basic::options::OptionKeys::rna::denovo::out::output_jump_res ]() );
	set_output_jump_o3p_to_o5p( opts[ basic::options::OptionKeys::rna::denovo::out::output_jump_o3p_to_o5p ]() );
	set_output_jump_chainbreak( opts[ basic::options::OptionKeys::rna::denovo::out::output_jump_chainbreak ]() );
	set_output_rotation_vector( opts[ basic::options::OptionKeys::rna::denovo::out::output_rotation_vector ]() );
	if ( opts[ basic::options::OptionKeys::rna::denovo::out::target_xyz ].user() ) {
		runtime_assert( opts[ basic::options::OptionKeys::rna::denovo::out::target_xyz ]().size() == 3 );
		set_output_jump_target_xyz( core::Vector(
			opts[ basic::options::OptionKeys::rna::denovo::out::target_xyz ]()[1],
			opts[ basic::options::OptionKeys::rna::denovo::out::target_xyz ]()[2],
			opts[ basic::options::OptionKeys::rna::denovo::out::target_xyz ]()[3]) );
	}
	set_output_jump_reference_RT_string( opts[ basic::options::OptionKeys::rna::denovo::out::output_jump_reference_RT ]() );
	set_save_jump_histogram( opts[ basic::options::OptionKeys::rna::denovo::out::save_jump_histogram ]() );
	set_jump_histogram_boxsize( opts[ basic::options::OptionKeys::rna::denovo::out::jump_histogram_boxsize ]() );
	set_jump_histogram_binwidth( opts[ basic::options::OptionKeys::rna::denovo::out::jump_histogram_binwidth ]() );
	set_jump_histogram_binwidth_rotvector( opts[ basic::options::OptionKeys::rna::denovo::out::jump_histogram_binwidth_rotvector ]() );
	set_output_histogram_file( opts[ basic::options::OptionKeys::rna::denovo::out::output_histogram_file ]() );
	if ( output_score_file_.size() == 0 && ( output_score_frequency_ != 0 || output_jump_res_.size() > 0 ) ) {
		output_score_file_ = opts[ out::file::silent ]();
		std::string::size_type pos = output_score_file_.find( ".out", 0 );
		std::string const new_prefix = ".SCORES.txt";
		if ( pos == std::string::npos ) utility_exit_with_message( "If you want to output a running score file, specify -output_score_file" );
		output_score_file_.replace( pos, new_prefix.length(), new_prefix );
	}
	if ( output_histogram_file_.size() == 0 && save_jump_histogram_ ) {
		output_histogram_file_ = option[ out::file::silent ]();
		std::string::size_type pos = output_histogram_file_.find( ".out", 0 );
		std::string const new_prefix = ".HISTOGRAM.bin.gz";
		if ( pos == std::string::npos ) utility_exit_with_message( "If you want to output a histogram, specify -output_histogram_file" );
		output_histogram_file_.replace( pos, new_prefix.length(), new_prefix );
	}

	set_allow_consecutive_bulges( opts[ basic::options::OptionKeys::rna::denovo::allow_consecutive_bulges ]() ) ;
	set_allowed_bulge_res( opts[ basic::options::OptionKeys::rna::denovo::allowed_bulge_res ]() ) ;

	if ( opts[ basic::options::OptionKeys::stepwise::rna::VDW_rep_screen_info ].user() ) {
		set_filter_vdw( true );
		VDW_rep_screen_info_ = opts[basic::options::OptionKeys::stepwise::rna::VDW_rep_screen_info]();
		set_vdw_rep_screen_include_sidechains( opts[ basic::options::OptionKeys::rna::denovo::VDW_rep_screen_include_sidechains ]() );
	}
	set_gradual_constraints( opts[ basic::options::OptionKeys::rna::denovo::gradual_constraints ]() );
	set_ramp_rnp_vdw( opts[ basic::options::OptionKeys::rna::denovo::ramp_rnp_vdw ]() );
	set_grid_vdw_weight( opts[ basic::options::OptionKeys::rna::denovo::grid_vdw_weight ]() );

	set_convert_protein_centroid( opts[ basic::options::OptionKeys::rna::denovo::convert_protein_CEN ]() );

	set_rna_protein_docking( opts[ basic::options::OptionKeys::rna::denovo::rna_protein_docking ]() );
	set_rna_protein_docking_legacy( opts[ basic::options::OptionKeys::rna::denovo::rna_protein_docking_legacy ]() );
	set_rna_protein_docking_freq( opts[ basic::options::OptionKeys::rna::denovo::rna_protein_docking_freq ]() );

	set_small_docking_moves( opts[ basic::options::OptionKeys::rna::denovo::small_docking_moves ]() );

	set_docking_move_size( opts[ basic::options::OptionKeys::rna::denovo::docking_move_size ]() );

	set_randomize_init_rnp( opts[ basic::options::OptionKeys::rna::denovo::randomize_init_rnp ]() );

	set_rnp_high_res_relax( opts[ basic::options::OptionKeys::rna::denovo::rnp_high_res_relax ]() );
	set_rnp_high_res_cycles( opts[ basic::options::OptionKeys::rna::denovo::rnp_high_res_cycles ]() );

	set_rnp_pack_first( opts[ basic::options::OptionKeys::rna::denovo::rnp_pack_first ]() );
	set_rnp_ramp_rep( opts[ basic::options::OptionKeys::rna::denovo::rnp_ramp_rep ]() );
	set_rnp_min_first( opts[ basic::options::OptionKeys::rna::denovo::rnp_min_first ]() );

	// dock into density is default false, can only be true if user specifies true AND a mapfile is provided
	if ( opts[ basic::options::OptionKeys::edensity::mapfile ].user() ) {
		set_dock_into_density_legacy( opts[ basic::options::OptionKeys::rna::denovo::dock_into_density_legacy ]() );
	}

	if ( rna_protein_docking() || dock_into_density() || opts[ basic::options::OptionKeys::rna::denovo::virtual_anchor ]().size() > 0 ) {
		set_docking( true );
	}

	set_new_fold_tree_initializer( opts[ basic::options::OptionKeys::rna::denovo::new_fold_tree_initializer ]() );

	if ( opts[ basic::options::OptionKeys::rna::denovo::initial_structures ].user() ) {
		set_initial_structures_provided( true );
	}

	set_simple_rmsd_cutoff_relax( opts[ basic::options::OptionKeys::rna::denovo::simple_relax ] );

	set_override_refine_pose_rounds( opts[ basic::options::OptionKeys::rna::denovo::override_refine_pose_rounds ]() );

	set_refine_native_get_good_FT( opts[ basic::options::OptionKeys::rna::denovo::refine_native_get_good_FT ]() );

	set_superimpose_over_all( opts[ stepwise::superimpose_over_all ]() );

	set_fixed_stems( opts[ basic::options::OptionKeys::rna::denovo::fixed_stems ]() );

	std::string const in_path = opts[ in::path::path ]()[1];
	if ( opts[ basic::options::OptionKeys::rna::denovo::vall_torsions ].user() ) {
		// check in database first
		std::string vall_torsions_file = basic::database::full_name("/sampling/rna/" + opts[ basic::options::OptionKeys::rna::denovo::vall_torsions ]() );
		if ( !utility::file::file_exists( vall_torsions_file ) && !utility::file::file_exists( vall_torsions_file + ".gz" ) ) vall_torsions_file = in_path + opts[ basic::options::OptionKeys::rna::denovo::vall_torsions ]();
		set_vall_torsions_file( vall_torsions_file );
	}
	if ( opts[ basic::options::OptionKeys::rna::denovo::use_1jj2_torsions ]() ) set_vall_torsions_file( basic::database::full_name("sampling/rna/1jj2.torsions") );
	if ( opts[ basic::options::OptionKeys::rna::denovo::jump_library_file ].user() ) set_jump_library_file( in_path + opts[ basic::options::OptionKeys::rna::denovo::jump_library_file] );
	if ( opts[ basic::options::OptionKeys::rna::denovo::params_file ].user() ) set_rna_params_file( in_path + opts[ basic::options::OptionKeys::rna::denovo::params_file ] );

	if ( opts[ basic::options::OptionKeys::rna::denovo::rna_lores_linear_chainbreak_weight ].user() ) set_linear_chainbreak_weight( opts[ basic::options::OptionKeys::rna::denovo::rna_lores_linear_chainbreak_weight ]() );
	if ( opts[ basic::options::OptionKeys::rna::denovo::rna_lores_chainbreak_weight ].user() ) set_chainbreak_weight( opts[ basic::options::OptionKeys::rna::denovo::rna_lores_chainbreak_weight ]() );

	if ( filter_lores_base_pairs_early_ ) set_filter_lores_base_pairs( true );

	if ( opts[ basic::options::OptionKeys::rna::denovo::no_filters ]() ) {
		set_autofilter( false );
		set_filter_chain_closure( false );
		set_filter_chain_closure_distance( false );
		set_filter_chain_closure_halfway( false );
		set_filter_lores_base_pairs( false );
		set_filter_lores_base_pairs_early( false );
	}

	// AMW: this doesn't correspond to an actual option. But, this way,
	// we will only do this song and dance if we have an align_pdb specified.
	disallow_realign_ = !opts[ stepwise::align_pdb ].user();
	rmsd_screen_ = opts[ stepwise::rmsd_screen ]();
	if ( opts[ stepwise::align_pdb ].user() ) {
		align_pdb_ = opts[ stepwise::align_pdb ]();
	} else {
		align_pdb_ = "";
	}

	if ( opts[ basic::options::OptionKeys::rna::denovo::symm_hack_arity ].user() ) set_symm_hack_arity( opts[ basic::options::OptionKeys::rna::denovo::symm_hack_arity ] );
}

void
RNA_FragmentMonteCarloOptions::list_options_read( utility::options::OptionKeyList & opts ) {
	using namespace basic::options;
	RNA_MinimizerOptions::list_options_read( opts ); // includes RNA_BasicOptions

	opts + OptionKeys::rna::denovo::cycles
		+ OptionKeys::rna::denovo::rounds
		+ OptionKeys::rna::denovo::minimize_rna
		+ OptionKeys::rna::denovo::relax_rna
		+ OptionKeys::rna::denovo::allow_bulge
		+ OptionKeys::rna::denovo::temperature
		+ OptionKeys::rna::denovo::ignore_secstruct
		+ OptionKeys::rna::denovo::jump_change_frequency
		+ OptionKeys::rna::denovo::close_loops
		+ OptionKeys::rna::denovo::close_loops_after_each_move
		+ OptionKeys::rna::denovo::staged_constraints
		+ OptionKeys::rna::denovo::filter_lores_base_pairs
		+ OptionKeys::rna::denovo::filter_lores_base_pairs_early
		+ OptionKeys::rna::denovo::suppress_bp_constraint
		+ OptionKeys::rna::denovo::filter_chain_closure
		+ OptionKeys::rna::denovo::filter_chain_closure_distance
		+ OptionKeys::rna::denovo::filter_chain_closure_halfway
		+ OptionKeys::rna::vary_geometry
		+ OptionKeys::rna::denovo::root_at_first_rigid_body
		+ OptionKeys::rna::denovo::autofilter
		+ OptionKeys::rna::denovo::bps_moves
		+ OptionKeys::rna::denovo::disallow_bps_at_extra_min_res
		+ OptionKeys::rna::denovo::allow_fragment_moves_in_bps
		+ OptionKeys::rna::denovo::frag_size
		+ OptionKeys::rna::denovo::out::output_score_frequency
		+ OptionKeys::rna::denovo::out::output_score_file
		+ OptionKeys::rna::denovo::out::output_jump_res
		+ OptionKeys::rna::denovo::out::output_jump_o3p_to_o5p
		+ OptionKeys::rna::denovo::out::output_jump_chainbreak
		+ OptionKeys::rna::denovo::out::output_rotation_vector
		+ OptionKeys::rna::denovo::out::target_xyz
		+ OptionKeys::rna::denovo::out::output_jump_reference_RT
		+ OptionKeys::rna::denovo::out::save_jump_histogram
		+ OptionKeys::rna::denovo::out::jump_histogram_boxsize
		+ OptionKeys::rna::denovo::out::jump_histogram_binwidth
		+ OptionKeys::rna::denovo::out::jump_histogram_binwidth_rotvector
		+ OptionKeys::rna::denovo::out::output_histogram_file
		+ OptionKeys::out::file::silent
		+ OptionKeys::rna::denovo::allow_consecutive_bulges
		+ OptionKeys::rna::denovo::allowed_bulge_res
		+ OptionKeys::stepwise::rna::VDW_rep_screen_info
		+ OptionKeys::rna::denovo::VDW_rep_screen_include_sidechains
		+ OptionKeys::rna::denovo::gradual_constraints
		+ OptionKeys::rna::denovo::grid_vdw_weight
		+ OptionKeys::rna::denovo::convert_protein_CEN
		+ OptionKeys::rna::denovo::rna_protein_docking
		+ OptionKeys::rna::denovo::rna_protein_docking_freq
		+ OptionKeys::rna::denovo::simple_relax
		+ OptionKeys::stepwise::superimpose_over_all
		+ OptionKeys::rna::denovo::fixed_stems
		+ OptionKeys::in::path::path
		+ OptionKeys::rna::denovo::vall_torsions
		+ OptionKeys::rna::denovo::use_1jj2_torsions
		+ OptionKeys::rna::denovo::jump_library_file
		+ OptionKeys::rna::denovo::params_file
		+ OptionKeys::rna::denovo::rna_lores_linear_chainbreak_weight
		+ OptionKeys::rna::denovo::rna_lores_chainbreak_weight
		+ OptionKeys::rna::denovo::no_filters
		+ OptionKeys::stepwise::align_pdb
		+ OptionKeys::stepwise::rmsd_screen
		+ OptionKeys::rna::denovo::dock_each_chunk
		+ OptionKeys::rna::denovo::dock_each_chunk_per_chain
		+ OptionKeys::rna::denovo::center_jumps_in_single_stranded
		+ OptionKeys::rna::denovo::ramp_rnp_vdw
		+ OptionKeys::rna::denovo::rna_protein_docking_freq
		+ OptionKeys::rna::denovo::small_docking_moves
		+ OptionKeys::rna::denovo::docking_move_size
		+ OptionKeys::rna::denovo::randomize_init_rnp
		+ OptionKeys::rna::denovo::rnp_high_res_relax
		+ OptionKeys::rna::denovo::rnp_high_res_cycles
		+ OptionKeys::rna::denovo::rnp_pack_first
		+ OptionKeys::rna::denovo::rnp_ramp_rep
		+ OptionKeys::rna::denovo::rnp_min_first
		+ OptionKeys::edensity::mapfile
		+ OptionKeys::rna::denovo::dock_into_density_legacy
		+ OptionKeys::rna::denovo::virtual_anchor
		+ OptionKeys::rna::denovo::new_fold_tree_initializer
		+ OptionKeys::rna::denovo::initial_structures
		+ OptionKeys::rna::denovo::override_refine_pose_rounds
		+ OptionKeys::rna::denovo::refine_native_get_good_FT
		+ OptionKeys::rna::denovo::symm_hack_arity;
}


//////////////////////////////////////////////////////////////////////////////////////////
void
RNA_FragmentMonteCarloOptions::initialize_for_farna_optimizer( Size const & cycles /* = 0 */ )
{
	initialize_from_command_line();
	if ( cycles > 0 ) set_monte_carlo_cycles( cycles );
	set_verbose( false );

	// following are different default values for RNA_DeNovoOptimizer, compared to rna_denovo defaults.
	if ( !option[ basic::options::OptionKeys::rna::denovo::bps_moves ].user() ) set_bps_moves( true );
	if ( !option[ basic::options::OptionKeys::rna::denovo::filter_lores_base_pairs ].user() ) set_filter_lores_base_pairs( false );
	if ( !option[ basic::options::OptionKeys::rna::denovo::filter_chain_closure ].user() ) set_filter_chain_closure( false );
}


} //options
} //denovo
} //rna
} //protocols
