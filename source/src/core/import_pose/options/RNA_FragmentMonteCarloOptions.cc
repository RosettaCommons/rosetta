// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/import_pose/options/RNA_FragmentMonteCarloOptions.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/import_pose/options/RNA_FragmentMonteCarloOptions.hh>

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
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.rna.denovo.options.RNA_FragmentMonteCarloOptions" );

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace utility::tag;

namespace core {
namespace import_pose {
namespace options {

//Constructor
RNA_FragmentMonteCarloOptions::RNA_FragmentMonteCarloOptions()
{}

/// @brief clone the options
RNA_FragmentMonteCarloOptionsOP
RNA_FragmentMonteCarloOptions::clone() const
{
	return utility::pointer::make_shared< RNA_FragmentMonteCarloOptions >( *this );
}

///////////////////////////////////////////////////////////////////
void
RNA_FragmentMonteCarloOptions::initialize_from_command_line() {
	RNA_MinimizerOptions::initialize_from_options( basic::options::option ); // includes RNA_BasicOptions
	initialize_from_options( basic::options::option );
}

void
RNA_FragmentMonteCarloOptions::initialize_from_options( utility::options::OptionCollection const & opts ) {
	RNA_MinimizerOptions::initialize_from_options( opts ); // includes RNA_BasicOptions

	if ( opts[ OptionKeys::rna::denovo::cycles ].user() ) {
		set_monte_carlo_cycles( opts[ OptionKeys::rna::denovo::cycles ]() );
		set_user_defined_cycles( true );
	}
	set_rounds( opts[ OptionKeys::rna::denovo::rounds ]() );
	set_minimize_structure( opts[ OptionKeys::rna::denovo::minimize_rna ]() );
	set_relax_structure( opts[ OptionKeys::rna::denovo::relax_rna ]() );
	set_allow_bulge( opts[ OptionKeys::rna::denovo::allow_bulge ]() );
	set_temperature( opts[ OptionKeys::rna::denovo::temperature ] );
	set_ignore_secstruct( opts[ OptionKeys::rna::denovo::ignore_secstruct ] );
	set_jump_change_frequency( opts[ OptionKeys::rna::denovo::jump_change_frequency ] );
	set_close_loops( opts[ OptionKeys::rna::denovo::close_loops] );
	set_close_loops_after_each_move( opts[ OptionKeys::rna::denovo::close_loops_after_each_move ] );
	set_staged_constraints( opts[ OptionKeys::rna::denovo::staged_constraints ] ) ;
	set_filter_lores_base_pairs( opts[ OptionKeys::rna::denovo::filter_lores_base_pairs] );
	set_filter_lores_base_pairs_early( opts[ OptionKeys::rna::denovo::filter_lores_base_pairs_early] );
	set_suppress_bp_constraint( opts[ OptionKeys::rna::denovo::suppress_bp_constraint ]() );
	set_filter_chain_closure( opts[ OptionKeys::rna::denovo::filter_chain_closure ]() );
	set_filter_chain_closure_distance( opts[ OptionKeys::rna::denovo::filter_chain_closure_distance ]() );
	set_filter_chain_closure_halfway( opts[ OptionKeys::rna::denovo::filter_chain_closure_halfway ]() );
	// Shouldn't set something from RNA_MinimizerOptions here -- then you have to keep
	// both definitions coordinated!!
	//set_vary_bond_geometry( opts[ OptionKeys::rna::vary_geometry ] );
	set_root_at_first_rigid_body( opts[ OptionKeys::rna::denovo::root_at_first_rigid_body ] );
	set_dock_each_chunk( opts[ OptionKeys::rna::denovo::dock_each_chunk ]() );
	set_dock_each_chunk_per_chain( opts[ OptionKeys::rna::denovo::dock_each_chunk_per_chain ]() );
	set_dock_chunks( opts[ OptionKeys::rna::denovo::dock_chunks ]() );
	set_center_jumps_in_single_stranded( opts[ OptionKeys::rna::denovo::center_jumps_in_single_stranded ]() );
	set_autofilter( opts[ OptionKeys::rna::denovo::autofilter ] );
	set_bps_moves( opts[ OptionKeys::rna::denovo::bps_moves ] );
	set_disallow_bps_at_extra_min_res( opts[ OptionKeys::rna::denovo::disallow_bps_at_extra_min_res ] );
	set_allow_fragment_moves_in_bps( opts[ OptionKeys::rna::denovo::allow_fragment_moves_in_bps ] );
	set_frag_size( opts[ OptionKeys::rna::denovo::frag_size ]() );
	set_use_chem_shift_data( opts[ OptionKeys::score::rna::rna_chemical_shift_exp_data].user() );

	set_output_score_frequency( opts[ OptionKeys::rna::denovo::out::output_score_frequency ]() );
	set_output_lores_silent_file( opts[ OptionKeys::rna::denovo::out::output_lores_silent_file ] );
	set_output_score_file( opts[ OptionKeys::rna::denovo::out::output_score_file ]() );
	set_output_jump_res( opts[ OptionKeys::rna::denovo::out::output_jump_res ]() );
	set_output_jump_o3p_to_o5p( opts[ OptionKeys::rna::denovo::out::output_jump_o3p_to_o5p ]() );
	set_output_jump_chainbreak( opts[ OptionKeys::rna::denovo::out::output_jump_chainbreak ]() );
	set_output_rotation_vector( opts[ OptionKeys::rna::denovo::out::output_rotation_vector ]() );
	if ( opts[ OptionKeys::rna::denovo::out::target_xyz ].user() ) {
		runtime_assert( opts[ OptionKeys::rna::denovo::out::target_xyz ]().size() == 3 );
		set_output_jump_target_xyz( core::Vector(
			opts[ OptionKeys::rna::denovo::out::target_xyz ]()[1],
			opts[ OptionKeys::rna::denovo::out::target_xyz ]()[2],
			opts[ OptionKeys::rna::denovo::out::target_xyz ]()[3]) );
	}
	set_output_jump_reference_RT_string( opts[ OptionKeys::rna::denovo::out::output_jump_reference_RT ]() );
	set_save_jump_histogram( opts[ OptionKeys::rna::denovo::out::save_jump_histogram ]() );
	set_jump_histogram_boxsize( opts[ OptionKeys::rna::denovo::out::jump_histogram_boxsize ]() );
	set_jump_histogram_binwidth( opts[ OptionKeys::rna::denovo::out::jump_histogram_binwidth ]() );
	set_jump_histogram_binwidth_rotvector( opts[ OptionKeys::rna::denovo::out::jump_histogram_binwidth_rotvector ]() );
	set_output_histogram_file( opts[ OptionKeys::rna::denovo::out::output_histogram_file ]() );
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

	set_allow_consecutive_bulges( opts[ OptionKeys::rna::denovo::allow_consecutive_bulges ]() ) ;
	set_allowed_bulge_res( opts[ OptionKeys::rna::denovo::allowed_bulge_res ]() ) ;

	set_output_filters(  opts[ OptionKeys::rna::denovo::output_filters ] );
	if ( opts[ OptionKeys::stepwise::rna::VDW_rep_screen_info ].user() ) {
		set_filter_vdw( true );
		VDW_rep_screen_info_ = opts[OptionKeys::stepwise::rna::VDW_rep_screen_info]();
		set_vdw_rep_screen_include_sidechains( opts[ OptionKeys::rna::denovo::VDW_rep_screen_include_sidechains ]() );
	}
	set_gradual_constraints( opts[ OptionKeys::rna::denovo::gradual_constraints ]() );
	set_ramp_rnp_vdw( opts[ OptionKeys::rna::denovo::ramp_rnp_vdw ]() );
	set_grid_vdw_weight( opts[ OptionKeys::rna::denovo::grid_vdw_weight ]() );

	set_convert_protein_centroid( opts[ OptionKeys::rna::denovo::convert_protein_CEN ]() );

	set_rna_protein_docking( opts[ OptionKeys::rna::denovo::rna_protein_docking ]() );
	set_rna_protein_docking_legacy( opts[ OptionKeys::rna::denovo::rna_protein_docking_legacy ]() );
	set_rna_protein_docking_freq( opts[ OptionKeys::rna::denovo::rna_protein_docking_freq ]() );

	set_small_docking_moves( opts[ OptionKeys::rna::denovo::small_docking_moves ]() );

	set_docking_move_size( opts[ OptionKeys::rna::denovo::docking_move_size ]() );

	set_randomize_init_rnp( opts[ OptionKeys::rna::denovo::randomize_init_rnp ]() );

	set_rnp_high_res_relax( opts[ OptionKeys::rna::denovo::rnp_high_res_relax ]() );
	set_rnp_high_res_cycles( opts[ OptionKeys::rna::denovo::rnp_high_res_cycles ]() );

	set_rnp_pack_first( opts[ OptionKeys::rna::denovo::rnp_pack_first ]() );
	set_rnp_ramp_rep( opts[ OptionKeys::rna::denovo::rnp_ramp_rep ]() );
	set_rnp_min_first( opts[ OptionKeys::rna::denovo::rnp_min_first ]() );

	// dock into density is default false, can only be true if user specifies true AND a mapfile is provided
	if ( opts[ OptionKeys::edensity::mapfile ].user() ) {
		set_dock_into_density_legacy( opts[ OptionKeys::rna::denovo::dock_into_density_legacy ]() );
	}

	if ( rna_protein_docking() || dock_into_density() || opts[ OptionKeys::rna::denovo::virtual_anchor ]().size() > 0 ) {
		set_docking( true );
	}

	set_new_fold_tree_initializer( opts[ OptionKeys::rna::denovo::new_fold_tree_initializer ]() );

	if ( opts[ OptionKeys::rna::denovo::initial_structures ].user() ) {
		set_initial_structures_provided( true );
	}

	set_simple_rmsd_cutoff_relax( opts[ OptionKeys::rna::denovo::simple_relax ] );

	set_override_refine_pose_rounds( opts[ OptionKeys::rna::denovo::override_refine_pose_rounds ]() );

	set_refine_native_get_good_FT( opts[ OptionKeys::rna::denovo::refine_native_get_good_FT ]() );

	set_superimpose_over_all( opts[ stepwise::superimpose_over_all ]() );

	set_fixed_stems( opts[ OptionKeys::rna::denovo::fixed_stems ]() );

	set_ft_close_chains( opts[ OptionKeys::rna::denovo::ft_close_chains ]() );

	std::string const in_path = opts[ in::path::path ]()[1];

	// vall torsions: check in database first
	if ( opts[ OptionKeys::rna::denovo::use_1jj2_torsions ]() ) {
		set_vall_torsions_file( basic::database::full_name("sampling/rna/1jj2.torsions") );
	} else {
		std::string vall_torsions_file = basic::database::full_name("sampling/rna/" + opts[ OptionKeys::rna::denovo::vall_torsions ]() );
		if ( !utility::file::file_exists( vall_torsions_file ) && !utility::file::file_exists( vall_torsions_file + ".gz" ) ) vall_torsions_file = in_path + opts[ OptionKeys::rna::denovo::vall_torsions ]();
		set_vall_torsions_file( vall_torsions_file );
	}

	std::string jump_library_file = basic::database::full_name("sampling/rna/" + opts[ OptionKeys::rna::denovo::jump_library_file ]() );
	if ( !utility::file::file_exists( jump_library_file ) && !utility::file::file_exists( jump_library_file+ ".gz" ) ) jump_library_file = in_path + opts[ OptionKeys::rna::denovo::jump_library_file ]();
	set_jump_library_file( jump_library_file );

	if ( opts[ OptionKeys::rna::denovo::rna_lores_linear_chainbreak_weight ].user() ) set_linear_chainbreak_weight( opts[ OptionKeys::rna::denovo::rna_lores_linear_chainbreak_weight ]() );
	if ( opts[ OptionKeys::rna::denovo::rna_lores_chainbreak_weight ].user() ) set_chainbreak_weight( opts[ OptionKeys::rna::denovo::rna_lores_chainbreak_weight ]() );

	if ( filter_lores_base_pairs_early_ ) set_filter_lores_base_pairs( true );

	if ( opts[ OptionKeys::rna::denovo::no_filters ]() ) {
		set_autofilter( false );
		set_filter_chain_closure( false );
		set_filter_chain_closure_distance( false ); // this is nuts b/c it is not a boolean?
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

	if ( opts[ OptionKeys::rna::denovo::symm_hack_arity ].user() ) set_symm_hack_arity( opts[ OptionKeys::rna::denovo::symm_hack_arity ]() );

	if ( opts[ OptionKeys::rna::denovo::exhaustive_fragment_insertion ].user() ) set_exhaustive_fragment_insertion( opts[ OptionKeys::rna::denovo::exhaustive_fragment_insertion ] );

	save_times_ = opts[ OptionKeys::out::save_times ]();
}


void
RNA_FragmentMonteCarloOptions::initialize_from_tag( utility::tag::TagCOP const & tag ) {
	RNA_MinimizerOptions::initialize_from_tag( tag ); // includes RNA_BasicOptions

	if ( tag->hasOption( "cycles" ) ) {
		set_monte_carlo_cycles( tag->getOption< core::Size >( "cycles", monte_carlo_cycles_ ) );
		set_user_defined_cycles( true );
	}

	if ( tag->hasOption( "rounds" ) ) {
		set_rounds( tag->getOption< core::Size >( "rounds", rounds_ ) );
	}

	if ( tag->hasOption( "minimize_rna" ) ) {
		set_minimize_structure( tag->getOption< bool >( "minimize_rna", minimize_structure_ ) );
	}

	if ( tag->hasOption( "relax_rna" ) ) {
		set_relax_structure( tag->getOption< bool >( "relax_rna", relax_structure_ ) );
	}

	if ( tag->hasOption( "allow_bulge" ) ) {
		set_allow_bulge( tag->getOption< bool >( "allow_bulge", allow_bulge_ ) );
	}

	if ( tag->hasOption( "temperature" ) ) {
		set_temperature( tag->getOption< core::Real >( "temperature", temperature_ ) );
	}

	if ( tag->hasOption( "ignore_secstruct" ) ) {
		set_ignore_secstruct( tag->getOption< bool >( "ignore_secstruct", ignore_secstruct_ ) );
	}

	if ( tag->hasOption( "jump_change_frequency" ) ) {
		set_jump_change_frequency( tag->getOption< core::Real >( "jump_change_frequency", jump_change_frequency_ ) );
	}

	if ( tag->hasOption( "close_loops" ) ) {
		set_close_loops( tag->getOption< bool >( "close_loops", close_loops_ ) );
	}

	if ( tag->hasOption( "close_loops_after_each_move" ) ) {
		set_close_loops_after_each_move( tag->getOption< bool >( "close_loops_after_each_move", close_loops_after_each_move_ ) );
	}

	if ( tag->hasOption( "staged_constraints" ) ) {
		set_staged_constraints( tag->getOption< bool >( "staged_constraints", staged_constraints_ ) );
	}

	if ( tag->hasOption( "filter_lores_base_pairs" ) ) {
		set_filter_lores_base_pairs( tag->getOption< bool >( "filter_lores_base_pairs", filter_lores_base_pairs_ ) );
	}

	if ( tag->hasOption( "filter_lores_base_pairs_early" ) ) {
		set_filter_lores_base_pairs_early( tag->getOption< bool >( "filter_lores_base_pairs_early", filter_lores_base_pairs_early_ ) );
	}

	if ( tag->hasOption( "suppress_bp_constraint" ) ) {
		set_suppress_bp_constraint( tag->getOption< core::Real >( "suppress_bp_constraint", suppress_bp_constraint_ ) );
	}

	if ( tag->hasOption( "filter_chain_closure" ) ) {
		set_filter_chain_closure( tag->getOption< bool >( "filter_chain_closure", filter_chain_closure_ ) );
	}

	if ( tag->hasOption( "filter_chain_closure_distance" ) ) {
		set_filter_chain_closure_distance( tag->getOption< core::Real >( "filter_chain_closure_distance", filter_chain_closure_distance_ ) );
	}

	if ( tag->hasOption( "filter_chain_closure_halfway" ) ) {
		set_filter_chain_closure_halfway( tag->getOption< bool >( "filter_chain_closure_halfway", filter_chain_closure_halfway_ ) );
	}

	if ( tag->hasOption( "root_at_first_rigid_body" ) ) {
		set_root_at_first_rigid_body( tag->getOption< bool >( "root_at_first_rigid_body", root_at_first_rigid_body_ ) );
	}

	if ( tag->hasOption( "dock_each_chunk" ) ) {
		set_dock_each_chunk( tag->getOption< bool >( "dock_each_chunk", dock_each_chunk_ ) );
	}

	if ( tag->hasOption( "dock_each_chunk_per_chain" ) ) {
		set_dock_each_chunk_per_chain( tag->getOption< bool >( "dock_each_chunk_per_chain", dock_each_chunk_per_chain_ ) );
	}

	if ( tag->hasOption( "dock_chunks" ) ) {
		auto dock_chunks_string = tag->getOption< std::string >( "dock_chunks", "" );
		set_dock_chunks( utility::string_split( dock_chunks_string ) );
	}

	if ( tag->hasOption( "center_jumps_in_single_stranded" ) ) {
		set_center_jumps_in_single_stranded( tag->getOption< bool >( "center_jumps_in_single_stranded", center_jumps_in_single_stranded_ ) );
	}

	if ( tag->hasOption( "autofilter" ) ) {
		set_autofilter( tag->getOption< bool >( "autofilter", autofilter_ ) );
	}

	if ( tag->hasOption( "bps_moves" ) ) {
		set_bps_moves( tag->getOption< bool >( "bps_moves", bps_moves_ ) );
	}

	if ( tag->hasOption( "disallow_bps_at_extra_min_res" ) ) {
		set_disallow_bps_at_extra_min_res( tag->getOption< bool >( "disallow_bps_at_extra_min_res", disallow_bps_at_extra_min_res_ ) );
	}

	if ( tag->hasOption( "allow_fragment_moves_in_bps" ) ) {
		set_allow_fragment_moves_in_bps( tag->getOption< bool >( "allow_fragment_moves_in_bps", allow_fragment_moves_in_bps_ ) );
	}

	if ( tag->hasOption( "frag_size" ) ) {
		set_frag_size( tag->getOption< core::Size >( "frag_size", frag_size_ ) );
	}

	// This is actually more complicated. This HAS to be provided on the command line.
	// Let's not touch this for now, but it's an AMW TODO.
	set_use_chem_shift_data( option[ OptionKeys::score::rna::rna_chemical_shift_exp_data].user() );
	/*
	if ( tag->hasOption( "rna_chemical_shift_exp_data" ) {
	set_use_chem_shift_data( true ); tag->getOption< bool >( "rna_chemical_shift_exp_data", use_chem_shift_data_ ) );
	}
	*/

	if ( tag->hasOption( "output_score_frequency" ) ) {
		set_output_score_frequency( tag->getOption< core::Real >( "output_score_frequency", output_score_frequency_ ) );
	}

	if ( tag->hasOption( "output_lores_silent_file" ) ) {
		set_output_lores_silent_file( tag->getOption< bool >( "output_lores_silent_file", output_lores_silent_file_ ) );
	}

	if ( tag->hasOption( "output_score_file" ) ) {
		set_output_score_file( tag->getOption< std::string >( "output_score_file", output_score_file_ ) );
	}

	if ( tag->hasOption( "output_jump_res" ) ) {
		std::istringstream ss;
		ss.str( tag->getOption< std::string >( "output_jump_res", "" ) );
		utility::vector1< core::Size > vec;
		while ( ss ) {
			core::Size n;
			ss >> n;
			vec.push_back( n );
		}
		set_output_jump_res( vec );
	}

	if ( tag->hasOption( "output_jump_o3p_to_o5p" ) ) {
		set_output_jump_o3p_to_o5p( tag->getOption< bool >( "output_jump_o3p_to_o5p", output_jump_o3p_to_o5p_ ) );
	}

	if ( tag->hasOption( "output_jump_chainbreak" ) ) {
		set_output_jump_chainbreak( tag->getOption< bool >( "output_jump_chainbreak", output_jump_chainbreak_ ) );
	}

	if ( tag->hasOption( "output_rotation_vector" ) ) {
		set_output_rotation_vector( tag->getOption< bool >( "output_rotation_vector", output_rotation_vector_ ) );
	}

	if ( tag->hasOption( "target_xyz" ) ) {
		// three reals, let's say space separated.
		std::istringstream ss;
		core::Real a, b, c;
		ss.str( tag->getOption< std::string >( "target_xyz", "" ) );
		ss >> a >> b >> c;
		set_output_jump_target_xyz( core::Vector( a, b, c ) );
	}

	if ( tag->hasOption( "output_jump_reference_RT" ) ) {
		set_output_jump_reference_RT_string( tag->getOption< std::string >( "output_jump_reference_RT", output_jump_reference_RT_string_ ) );
	}

	if ( tag->hasOption( "save_jump_histogram" ) ) {
		set_save_jump_histogram( tag->getOption< bool >( "save_jump_histogram", save_jump_histogram_ ) );
	}

	if ( tag->hasOption( "jump_histogram_boxsize" ) ) {
		set_jump_histogram_boxsize( tag->getOption< core::Real >( "jump_histogram_boxsize", jump_histogram_boxsize_ ) );
	}

	if ( tag->hasOption( "jump_histogram_binwidth" ) ) {
		set_jump_histogram_binwidth( tag->getOption< core::Real >( "jump_histogram_binwidth", jump_histogram_binwidth_ ) );
	}

	if ( tag->hasOption( "jump_histogram_binwidth_rotvector" ) ) {
		set_jump_histogram_binwidth_rotvector( tag->getOption< core::Real >( "jump_histogram_binwidth_rotvector", jump_histogram_binwidth_rotvector_ ) );
	}

	if ( tag->hasOption( "output_histogram_file" ) ) {
		set_output_histogram_file( tag->getOption< std::string >( "output_histogram_file", output_histogram_file_ ) );
	}

	if ( output_score_file_.size() == 0 && ( output_score_frequency_ != 0 || output_jump_res_.size() > 0 ) ) {
		utility_exit_with_message( "If you want to output a running score file, specify -output_score_file" );
	}
	if ( output_histogram_file_.size() == 0 && save_jump_histogram_ ) {
		utility_exit_with_message( "If you want to output a histogram, specify -output_histogram_file" );
	}

	if ( tag->hasOption( "allow_consecutive_bulges" ) ) {
		set_allow_consecutive_bulges( tag->getOption< bool >( "allow_consecutive_bulges", allow_consecutive_bulges_ ) );
	}

	if ( tag->hasOption( "allowed_bulge_res" ) ) {
		std::istringstream ss;
		ss.str( tag->getOption< std::string >( "allowed_bulge_res", "" ) );
		utility::vector1< core::Size > vec;
		while ( ss ) {
			core::Size n;
			ss >> n;
			vec.push_back( n );
		}
		set_allowed_bulge_res( vec );
	}

	if ( tag->hasOption( "output_filters" ) ) {
		set_output_filters( tag->getOption< bool >( "output_filters", output_filters_ ) );
	}

	if ( tag->hasOption( "VDW_rep_screen_info" ) ) {
		set_filter_vdw( true );
		set_VDW_rep_screen_info( utility::string_split( tag->getOption< std::string >( "VDW_rep_screen_info", "" ) ) );
		if ( tag->hasOption( "vdw_rep_screen_include_sidechains" ) ) {
			set_vdw_rep_screen_include_sidechains( tag->getOption< bool >( "vdw_rep_screen_include_sidechains", vdw_rep_screen_include_sidechains_ ) );
		}
	}

	if ( tag->hasOption( "gradual_constraints" ) ) {
		set_gradual_constraints( tag->getOption< bool >( "gradual_constraints", gradual_constraints_ ) );
	}

	if ( tag->hasOption( "ramp_rnp_vdw" ) ) {
		set_ramp_rnp_vdw( tag->getOption< bool >( "ramp_rnp_vdw", ramp_rnp_vdw_ ) );
	}

	if ( tag->hasOption( "grid_vdw_weight" ) ) {
		set_grid_vdw_weight( tag->getOption< core::Real >( "grid_vdw_weight", grid_vdw_weight_ ) );
	}

	if ( tag->hasOption( "convert_protein_centroid" ) ) {
		set_convert_protein_centroid( tag->getOption< bool >( "convert_protein_centroid", convert_protein_centroid_ ) );
	}

	if ( tag->hasOption( "rna_protein_docking" ) ) {
		set_rna_protein_docking( tag->getOption< bool >( "rna_protein_docking", rna_protein_docking_ ) );
	}
	if ( tag->hasOption( "rna_protein_docking_legacy" ) ) {
		set_rna_protein_docking_legacy( tag->getOption< bool >( "rna_protein_docking_legacy", rna_protein_docking_legacy_ ) );
	}
	if ( tag->hasOption( "rna_protein_docking_freq" ) ) {
		set_rna_protein_docking_freq( tag->getOption< core::Real >( "rna_protein_docking_freq", rna_protein_docking_freq_ ) );
	}

	if ( tag->hasOption( "small_docking_moves" ) ) {
		set_small_docking_moves( tag->getOption< bool >( "small_docking_moves", small_docking_moves_ ) );
	}
	if ( tag->hasOption( "docking_move_size" ) ) {
		set_docking_move_size( tag->getOption< core::Real >( "docking_move_size", docking_move_size_ ) );
	}

	if ( tag->hasOption( "randomize_init_rnp" ) ) {
		set_randomize_init_rnp( tag->getOption< bool >( "randomize_init_rnp", randomize_init_rnp_ ) );
	}
	if ( tag->hasOption( "rnp_high_res_relax" ) ) {
		set_rnp_high_res_relax( tag->getOption< bool >( "rnp_high_res_relax", rnp_high_res_relax_ ) );
	}
	if ( tag->hasOption( "rnp_high_res_cycles" ) ) {
		set_rnp_high_res_cycles( tag->getOption< core::Size >( "rnp_high_res_cycles", rnp_high_res_cycles_ ) );
	}

	if ( tag->hasOption( "rnp_pack_first" ) ) {
		set_rnp_pack_first( tag->getOption< bool >( "rnp_pack_first", rnp_pack_first_ ) );
	}
	if ( tag->hasOption( "rnp_ramp_rep" ) ) {
		set_rnp_ramp_rep( tag->getOption< bool >( "rnp_ramp_rep", rnp_ramp_rep_ ) );
	}
	if ( tag->hasOption( "rnp_min_first" ) ) {
		set_rnp_min_first( tag->getOption< bool >( "rnp_min_first", rnp_min_first_ ) );
	}

	// UGH WE NEED TO REALLY CHECK OPTIONS HERE AMW TODO #3627
	// dock into density is default false, can only be true if user specifies true AND a mapfile is provided
	if ( option[ OptionKeys::edensity::mapfile ].user() ) {
		if ( tag->hasOption( "dock_into_density_legacy" ) ) {
			set_dock_into_density_legacy( tag->getOption< bool >( "dock_into_density_legacy", dock_into_density_legacy_ ) );
		}
	}

	// Read in virtual anchor just to check... really could just check if it's specified.
	if ( rna_protein_docking() || dock_into_density() || tag->hasOption( "virtual_anchor" ) ) {
		set_docking( true );
	}

	if ( tag->hasOption( "new_fold_tree_initializer" ) ) {
		set_new_fold_tree_initializer( tag->getOption< bool >( "new_fold_tree_initializer", new_fold_tree_initializer_ ) );
	}

	if ( tag->hasOption( "initial_structures" ) ) {
		set_initial_structures_provided( true );
	}

	if ( tag->hasOption( "simple_relax" ) ) {
		set_simple_rmsd_cutoff_relax( tag->getOption< bool >( "simple_relax", simple_rmsd_cutoff_relax_ ) );
	}
	if ( tag->hasOption( "override_refine_pose_rounds" ) ) {
		set_override_refine_pose_rounds( tag->getOption< bool >( "override_refine_pose_rounds", override_refine_pose_rounds_ ) );
	}
	if ( tag->hasOption( "refine_native_get_good_FT" ) ) {
		set_refine_native_get_good_FT( tag->getOption< bool >( "refine_native_get_good_FT", refine_native_get_good_FT_ ) );
	}
	if ( tag->hasOption( "superimpose_over_all" ) ) {
		set_superimpose_over_all( tag->getOption< bool >( "superimpose_over_all", superimpose_over_all_ ) );
	}
	if ( tag->hasOption( "fixed_stems" ) ) {
		set_fixed_stems( tag->getOption< bool >( "fixed_stems", fixed_stems_ ) );
	}
	if ( tag->hasOption( "ft_close_chains" ) ) {
		set_ft_close_chains( tag->getOption< bool >( "ft_close_chains", ft_close_chains_ ) );
	}

	// AMW TODO: check
	if ( tag->hasOption( "use_1jj2_torsions" ) ) {
		set_vall_torsions_file( basic::database::full_name("sampling/rna/1jj2.torsions") );
	} else {
		std::string vall_torsions = basic::database::full_name("sampling/rna/"
			+ tag->getOption< std::string >( "vall_torsions", "RNA18_HUB_2.154_2.5.torsions" ) );
		if ( !utility::file::file_exists( vall_torsions ) && !utility::file::file_exists( vall_torsions+ ".gz" ) ) {
			vall_torsions = "./"
				+ tag->getOption< std::string >( "vall_torsions", "RNA18_HUB_2.154_2.5.torsions" );
		}
		set_vall_torsions_file( vall_torsions );
	}

	std::string jump_library_file = basic::database::full_name("sampling/rna/"
		+ tag->getOption< std::string >( "jump_library_file", "RNA18_HUB_2.154_2.5.jump" ) );
	if ( !utility::file::file_exists( jump_library_file ) && !utility::file::file_exists( jump_library_file+ ".gz" ) ) {
		jump_library_file = "./"
			+ tag->getOption< std::string >( "jump_library_file", "RNA18_HUB_2.154_2.5.jump" );
	}
	set_jump_library_file( jump_library_file );


	if ( tag->hasOption( "linear_chainbreak_weight" ) ) {
		set_linear_chainbreak_weight( tag->getOption< core::Real >( "linear_chainbreak_weight", linear_chainbreak_weight_ ) );
	}
	if ( tag->hasOption( "chainbreak_weight" ) ) {
		set_chainbreak_weight( tag->getOption< core::Real >( "chainbreak_weight", chainbreak_weight_ ) );
	}

	if ( filter_lores_base_pairs_early_ ) set_filter_lores_base_pairs( true );


	if ( tag->hasOption( "no_filters" ) ) {
		set_autofilter( false );
		set_filter_chain_closure( false );
		set_filter_chain_closure_distance( false );
		set_filter_chain_closure_halfway( false );
		set_filter_lores_base_pairs( false );
		set_filter_lores_base_pairs_early( false );
	}


	if ( tag->hasOption( "align_pdb" ) ) {
		set_disallow_realign( false );
		set_align_pdb( tag->getOption< std::string >( "align_pdb", align_pdb_ ) );
	}
	if ( tag->hasOption( "rmsd_screen" ) ) {
		set_rmsd_screen( tag->getOption< core::Real >( "rmsd_screen", rmsd_screen_ ) );
	}

	if ( tag->hasOption( "symm_hack_arity" ) ) {
		set_symm_hack_arity( tag->getOption< core::Size >( "symm_hack_arity", symm_hack_arity_ ) );
	}

	if ( tag->hasOption( "exhaustive_fragment_insertion" ) ) {
		set_exhaustive_fragment_insertion( tag->getOption< bool >( "exhaustive_fragment_insertion", exhaustive_fragment_insertion_ ) );
	}
	if ( tag->hasOption( "save_times" ) ) {
		set_save_times( tag->getOption< bool >( "save_times", save_times_ ) );
	}
}


void
RNA_FragmentMonteCarloOptions::list_attributes( AttributeList & attlist ) {

	attlist
		+ XMLSchemaAttribute::attribute_w_default( "cycles", xsct_non_negative_integer, "Number of Monte Carlo cycles", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "rounds", xsct_non_negative_integer, "Number of Monte Carlo rounds", "10" )
		+ XMLSchemaAttribute::required_attribute( "minimize_rna", xsct_rosetta_bool, "XRW TODO" )
		+ XMLSchemaAttribute::attribute_w_default( "relax_rna", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "allow_bulge", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "temperature", xsct_real, "XRW TODO", "2.0" )
		+ XMLSchemaAttribute::attribute_w_default( "ignore_secstruct", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "jump_change_frequency", xsct_real, "XRW TODO", "0.1" )
		+ XMLSchemaAttribute::attribute_w_default( "close_loops", xsct_rosetta_bool, "XRW TODO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "close_loops_after_each_move", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "staged_constraints", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "filter_lores_base_pairs", xsct_rosetta_bool, "XRW TODO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "filter_lores_base_pairs_early", xsct_rosetta_bool, "XRW TODO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "suppress_bp_constraint", xsct_real, "XRW TODO", "1.0" )
		+ XMLSchemaAttribute::attribute_w_default( "filter_chain_closure", xsct_rosetta_bool, "XRW TODO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "filter_chain_closure_distance", xsct_real, "XRW TODO", "6.0" )
		+ XMLSchemaAttribute::attribute_w_default( "filter_chain_closure_halfway", xsct_rosetta_bool, "XRW TODO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "root_at_first_rigid_body", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "dock_each_chunk", xsct_rosetta_bool, "XRW TODO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "dock_each_chunk_per_chain", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "dock_chunks", xs_string, "space separated PDB list", "" )
		+ XMLSchemaAttribute::attribute_w_default( "center_jumps_in_single_stranded", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "autofilter", xsct_rosetta_bool, "XRW TODO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "bps_moves", xsct_rosetta_bool, "XRW TODO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "disallow_bps_at_extra_min_res", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "allow_fragment_moves_in_bps", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "frag_size", xsct_non_negative_integer, "XRW TODO", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "output_score_frequency", xsct_real, "XRW TODO", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "output_lores_silent_file", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "output_score_file", xs_string, "XRW TODO", "" )
		+ XMLSchemaAttribute::attribute_w_default( "output_jump_res", xsct_nnegative_int_wsslist, "XRW TODO", "" )
		+ XMLSchemaAttribute::attribute_w_default( "output_jump_o3p_to_o5p", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "output_jump_chainbreak", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "output_rotation_vector", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "target_xyz", xs_string, "AMW TODO", "" )
		+ XMLSchemaAttribute::attribute_w_default( "output_jump_reference_RT", xs_string, "XRW TODO", "" )
		+ XMLSchemaAttribute::attribute_w_default( "save_jump_histogram", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "jump_histogram_boxsize", xsct_real, "XRW TODO", "40.0" )
		+ XMLSchemaAttribute::attribute_w_default( "jump_histogram_binwidth", xsct_real, "XRW TODO", "4.0" )
		+ XMLSchemaAttribute::attribute_w_default( "jump_histogram_binwidth_rotvector", xsct_real, "XRW TODO", "36.0" )
		+ XMLSchemaAttribute::attribute_w_default( "output_histogram_file", xs_string, "XRW TODO", "" )
		+ XMLSchemaAttribute::attribute_w_default( "allow_consecutive_bulges", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "allowed_bulge_res", xsct_nnegative_int_wsslist, "XRW TODO", "" )
		+ XMLSchemaAttribute::attribute_w_default( "output_filters", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "VDW_rep_screen_info", xs_string, "XRW TODO", "" )
		+ XMLSchemaAttribute::attribute_w_default( "vdw_rep_screen_include_sidechains", xsct_rosetta_bool, "XRW TODO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "gradual_constraints", xsct_rosetta_bool, "XRW TODO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "ramp_rnp_vdw", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "grid_vdw_weight", xsct_real, "XRW TODO", "1.0" )
		+ XMLSchemaAttribute::attribute_w_default( "convert_protein_centroid", xsct_rosetta_bool, "XRW TODO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "rna_protein_docking", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "rna_protein_docking_legacy", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "rna_protein_docking_freq", xsct_real, "XRW TODO", "0.4" )
		+ XMLSchemaAttribute::attribute_w_default( "small_docking_moves", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "docking_move_size", xsct_real, "XRW TODO", "1.0" )
		+ XMLSchemaAttribute::attribute_w_default( "randomize_init_rnp", xsct_rosetta_bool, "XRW TODO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "rnp_high_res_relax", xsct_rosetta_bool, "XRW TODO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "rnp_high_res_cycles", xsct_non_negative_integer, "XRW TODO", "10" )
		+ XMLSchemaAttribute::attribute_w_default( "rnp_pack_first", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "rnp_ramp_rep", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "rnp_min_first", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "dock_into_density_legacy", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "new_fold_tree_initializer", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "simple_relax", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "override_refine_pose_rounds", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "refine_native_get_good_FT", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "superimpose_over_all", xsct_rosetta_bool, "XRW TODO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "fixed_stems", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "ft_close_chains", xsct_rosetta_bool, "XRW TODO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "use_1jj2_torsions", xs_string, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "vall_torsions", xs_string, "XRW TODO", "RNA18_HUB_2.154_2.5.torsions" )
		+ XMLSchemaAttribute::attribute_w_default( "jump_library_file", xs_string, "XRW TODO", "RNA18_HUB_2.154_2.5.jump" )
		+ XMLSchemaAttribute::attribute_w_default( "linear_chainbreak_weight", xsct_real, "XRW TODO", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "chainbreak_weight", xsct_real, "XRW TODO", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "no_filters", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "rmsd_screen", xsct_real, "XRW TODO", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "symm_hack_arity", xsct_non_negative_integer, "XRW TODO", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "exhaustive_fragment_insertion", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "save_times", xsct_rosetta_bool, "XRW TODO", "false" );

	RNA_MinimizerOptions::list_attributes( attlist );
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
		+ OptionKeys::rna::denovo::out::output_lores_silent_file
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
		+ OptionKeys::rna::denovo::rna_protein_docking_legacy
		+ OptionKeys::rna::denovo::rna_protein_docking_freq
		+ OptionKeys::rna::denovo::simple_relax
		+ OptionKeys::stepwise::superimpose_over_all
		+ OptionKeys::rna::denovo::fixed_stems
		+ OptionKeys::rna::denovo::ft_close_chains
		+ OptionKeys::in::path::path
		+ OptionKeys::rna::denovo::vall_torsions
		+ OptionKeys::rna::denovo::use_1jj2_torsions
		+ OptionKeys::rna::denovo::jump_library_file
		+ OptionKeys::rna::denovo::rna_lores_linear_chainbreak_weight
		+ OptionKeys::rna::denovo::rna_lores_chainbreak_weight
		+ OptionKeys::rna::denovo::output_filters
		+ OptionKeys::rna::denovo::no_filters
		+ OptionKeys::stepwise::align_pdb
		+ OptionKeys::stepwise::rmsd_screen
		+ OptionKeys::rna::denovo::dock_each_chunk
		+ OptionKeys::rna::denovo::dock_each_chunk_per_chain
		+ OptionKeys::rna::denovo::dock_chunks
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
		+ OptionKeys::rna::denovo::helical_substructs
		+ OptionKeys::rna::denovo::override_refine_pose_rounds
		+ OptionKeys::rna::denovo::refine_native_get_good_FT
		+ OptionKeys::rna::denovo::symm_hack_arity
		+ OptionKeys::out::save_times
		+ OptionKeys::score::rna::rna_chemical_shift_exp_data
		+ OptionKeys::rna::denovo::exhaustive_fragment_insertion;
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
} //protocols
