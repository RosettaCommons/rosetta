// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/farna/options/RNA_FragmentMonteCarloOptions.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/farna/options/RNA_FragmentMonteCarloOptions.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh> // for superimpose_over_all
#include <basic/database/open.hh>
#include <utility/file/file_sys_util.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.farna.options.RNA_FragmentMonteCarloOptions" );

using namespace basic::options;
using namespace basic::options::OptionKeys;

namespace protocols {
namespace farna {
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
	chainbreak_weight_( -1.0 ), /* use farna/rna_lores.wts number unless user specified. -1.0 is never really used. */
	linear_chainbreak_weight_( -1.0 ),  /* use farna/rna_lores.wts number unless user specified. -1.0 is never really used. */
	close_loops_( true ),
	close_loops_in_last_round_( true ),
	close_loops_after_each_move_( false ),
	allow_bulge_( false ),
	allow_consecutive_bulges_( false ),
	jump_change_frequency_( 0.1 ),
	autofilter_( false ),
	autofilter_score_quantile_( 0.20 ),
	titrate_stack_bonus_( true ),
	root_at_first_rigid_body_( false ),
	suppress_bp_constraint_( 1.0 ),
	filter_lores_base_pairs_( false ),
	filter_lores_base_pairs_early_( false ),
	filter_chain_closure_( true ),
	filter_chain_closure_distance_( 6.0 ), /* in Angstroms. This is pretty loose!*/
	filter_chain_closure_halfway_( true ),
	staged_constraints_( false ),
	output_score_frequency_( 0 ),
	filter_vdw_( false ),
	vdw_rep_screen_include_sidechains_( false ),
	gradual_constraints_( true ),
	grid_vdw_weight_( 1.0 ),
	convert_protein_centroid_( true ),
	rna_protein_docking_( true ),
	rna_protein_docking_freq_( 10 ),
	simple_rmsd_cutoff_relax_( false ),
	refine_from_silent_( false ),
	refine_pose_( false ),
	bps_moves_( false ),
	disallow_bps_at_extra_min_res_( false ),
	allow_fragment_moves_in_bps_( false ),
	// following is odd, but note that core::scoring::rna::chemical_shift machinery also checks global options system.
	use_chem_shift_data_( option[ OptionKeys::score::rna_chemical_shift_exp_data].user() ),
	superimpose_over_all_( false ),
	fixed_stems_( false ),
	all_rna_fragments_file_( basic::database::full_name("sampling/rna/RICHARDSON_RNA09.torsions") ),
	rna_params_file_( "" ),
	jump_library_file_( basic::database::full_name("sampling/rna/1jj2_RNA_jump_library.dat" ) )
{}

//Destructor
RNA_FragmentMonteCarloOptions::~RNA_FragmentMonteCarloOptions()
{}

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

	RNA_MinimizerOptions::initialize_from_command_line(); // includes RNA_BasicOptions

	if ( option[ rna::farna::cycles ].user() ) {
		set_monte_carlo_cycles( option[ rna::farna::cycles ]() );
		set_user_defined_cycles( true );
	}
	set_rounds( option[ rna::farna::rounds ]() );
	minimize_structure_ = option[ rna::farna::minimize_rna ]();
	relax_structure_ = option[ rna::farna::relax_rna ]();
	allow_bulge_ = option[ rna::farna::allow_bulge ]();
	set_temperature( option[ rna::farna::temperature ] );
	set_ignore_secstruct( option[ rna::farna::ignore_secstruct ] );
	set_jump_change_frequency( option[ rna::farna::jump_change_frequency ] );
	set_close_loops( option[ rna::farna::close_loops] );
	set_close_loops_after_each_move( option[ rna::farna::close_loops_after_each_move ] );
	set_staged_constraints( option[ rna::farna::staged_constraints ] ) ;
	set_filter_lores_base_pairs(  option[ rna::farna::filter_lores_base_pairs] );
	set_filter_lores_base_pairs_early(  option[ rna::farna::filter_lores_base_pairs_early] );
	set_suppress_bp_constraint(  option[ rna::farna::suppress_bp_constraint]() );
	set_filter_chain_closure(  option[ rna::farna::filter_chain_closure ]() );
	set_filter_chain_closure_distance(  option[ rna::farna::filter_chain_closure_distance ]() );
	set_filter_chain_closure_halfway(  option[ rna::farna::filter_chain_closure_halfway ]() );
	set_vary_bond_geometry(  option[ rna::vary_geometry ] );
	set_root_at_first_rigid_body(  option[ rna::farna::root_at_first_rigid_body ] );
	set_autofilter( option[ rna::farna::autofilter ] );
	set_bps_moves( option[ rna::farna::bps_moves ] );
	set_disallow_bps_at_extra_min_res( option[ rna::farna::disallow_bps_at_extra_min_res ] );
	set_allow_fragment_moves_in_bps( option[ rna::farna::allow_fragment_moves_in_bps ] );

	set_output_score_frequency( option[ rna::farna::output_score_frequency ]() );
	set_output_score_file( option[ rna::farna::output_score_file ]() );
	if ( output_score_file_.size() == 0 && output_score_frequency_ != 0 )  {
		output_score_file_ = option[ out::file::silent ]();
		std::string::size_type pos = output_score_file_.find( ".out", 0 );
		std::string const new_prefix = ".SCORES.txt";
		if ( pos == std::string::npos ) utility_exit_with_message(  "If you want to output a running score file, specify -output_score_file" );
		output_score_file_.replace( pos, new_prefix.length(), new_prefix );
	}

	set_allow_consecutive_bulges( option[ rna::farna::allow_consecutive_bulges ]() ) ;
	set_allowed_bulge_res( option[ rna::farna::allowed_bulge_res ]() ) ;

	if ( option[ basic::options::OptionKeys::stepwise::rna::VDW_rep_screen_info ].user() ) {
		set_filter_vdw( true );
		set_vdw_rep_screen_include_sidechains( option[ rna::farna::VDW_rep_screen_include_sidechains ]() );
	}
	set_gradual_constraints( option[ rna::farna::gradual_constraints ]() );
	set_grid_vdw_weight( option[ rna::farna::grid_vdw_weight ]() );

	set_convert_protein_centroid( option[ rna::farna::convert_protein_CEN ]() ); // default true

	// default true, but of course only happens if structure contains both rna and protein
	set_rna_protein_docking( option[ rna::farna::rna_protein_docking ]() );
	set_rna_protein_docking_freq( option[ rna::farna::rna_protein_docking_freq ]() ); // default 10

	set_simple_rmsd_cutoff_relax( option[ rna::farna::simple_relax ] );

	set_superimpose_over_all( option[ stepwise::superimpose_over_all ]() );

	set_fixed_stems( option[ rna::farna::fixed_stems ]() );

	std::string const in_path = option[ in::path::path ]()[1];
	if ( option[ rna::farna::vall_torsions ].user() ) {
		// check in database first
		std::string vall_torsions_file = basic::database::full_name("/sampling/rna/" + option[ rna::farna::vall_torsions ]() );
		if ( !utility::file::file_exists( vall_torsions_file ) && !utility::file::file_exists( vall_torsions_file + ".gz" ) )  vall_torsions_file = in_path + option[ rna::farna::vall_torsions ]();
		set_vall_torsions_file( vall_torsions_file );
	}
	if ( option[ rna::farna::use_1jj2_torsions ]() ) set_vall_torsions_file( basic::database::full_name("sampling/rna/1jj2.torsions") );
	if ( option[ rna::farna::jump_library_file ].user() ) set_jump_library_file( in_path + option[ rna::farna::jump_library_file] );
	if ( option[ rna::farna::params_file ].user() ) set_rna_params_file( in_path + option[ rna::farna::params_file ] );

	if ( option[ rna::farna::rna_lores_linear_chainbreak_weight ].user() ) set_linear_chainbreak_weight( option[ rna::farna::rna_lores_linear_chainbreak_weight ]() );
	if ( option[ rna::farna::rna_lores_chainbreak_weight ].user() ) set_chainbreak_weight( option[ rna::farna::rna_lores_chainbreak_weight ]() );

	if ( filter_lores_base_pairs_early_ ) set_filter_lores_base_pairs( true );

	if ( option[ rna::farna::no_filters ]() ) {
		set_autofilter( false );
		set_filter_chain_closure( false );
		set_filter_chain_closure_distance( false );
		set_filter_chain_closure_halfway(  false );
		set_filter_lores_base_pairs(  false );
		set_filter_lores_base_pairs_early(  false );
	}
	
	// AMW: this doesn't correspond to an actual option. But, this way,
	// we will only do this song and dance if we have an align_pdb specified.
	disallow_realign_ = !option[ stepwise::align_pdb ].user();
	rmsd_screen_ = option[ stepwise::rmsd_screen ]();
	if ( option[ stepwise::align_pdb ].user() ) {
		align_pdb_ = option[ stepwise::align_pdb ]();
	} else {
		align_pdb_ = "";
	}
}


//////////////////////////////////////////////////////////////////////////////////////////
void
RNA_FragmentMonteCarloOptions::initialize_for_farna_optimizer( Size const & cycles /* = 0 */ )
{
	initialize_from_command_line();
	if ( cycles > 0 ) set_monte_carlo_cycles( cycles );
	set_verbose( false );

	// following are different default values for FARNA_Optimizer, compared to rna_denovo defaults.
	if ( !option[ rna::farna::bps_moves ].user() ) set_bps_moves( true );
	if ( !option[ rna::farna::filter_lores_base_pairs ].user() ) set_filter_lores_base_pairs( false );
	if ( !option[ rna::farna::filter_chain_closure ].user() ) set_filter_chain_closure( false );
}


} //options
} //farna
} //protocols
