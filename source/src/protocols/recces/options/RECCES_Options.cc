// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/options/RECCES_Options.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/recces/options/RECCES_Options.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/recces.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.recces.options.RECCES_Options" );

using namespace basic::options;
using namespace basic::options::OptionKeys;

namespace protocols {
namespace recces {
namespace options {

//Constructor
RECCES_Options::RECCES_Options():
	n_cycle_( 10000 ),
	n_dump_( 0 ),
	a_form_range_( 60.0 ),
	base_pair_rmsd_cutoff_( 0.0 ),
	base_pair_rotation_mag_( 0.0 ),
	base_pair_translation_mag_( 0.0 ),
	legacy_turner_mode_( false ),
	dump_pdb_( false ),
	save_scores_( false ),
	block_stack_( false ),
	dump_freq_( 500 ),
	dump_silent_( false ),
	out_torsions_( false ),
	setup_base_pair_constraints_( false ),
	seq1_( "" ),
	seq2_( "" ),
	infile_( "" ),
	out_prefix_( "" ),
	xyz_file_( "" ),
	histogram_min_( -100.05 ),
	histogram_max_( +800.05 ),
	histogram_spacing_( 0.1 ),
	angle_range_bb_( 20.0 ),
	angle_range_free_bb_( 180.0 ),
	angle_range_chi_( 20.0 ),
	angle_range_free_chi_( 180.0 ),
	chi_stdev_( 20.0 ),
	bb_stdev_( 1.0 ),
	standard_bb_stdev_( 1.0 ),
	thermal_sampler_mode_( false ),
	blank_score_terms_( false ),
	skip_last_accept_( false ),
	suppress_sampler_display_( false ),
	prefix_each_output_pdb_( false ),
	show_more_pose_scores_( false ),
	output_simple_text_files_( false ),
	accept_no_op_moves_( false ),
	sample_jump_( false )
{}

//Destructor
RECCES_Options::~RECCES_Options()
{}

/// @brief copy constructor
RECCES_Options::RECCES_Options( RECCES_Options const & src ) :
	ResourceOptions( src )
{
	*this = src;
}

/// @brief clone the options
RECCES_OptionsOP
RECCES_Options::clone() const
{
	return RECCES_OptionsOP( new RECCES_Options( *this ) );
}

///////////////////////////////////////////////////////////////////
void
RECCES_Options::initialize_from_command_line() {

	set_n_cycle( option[ OptionKeys::recces::n_cycle ]() );
	set_n_dump( option[ OptionKeys::recces::n_intermediate_dump ]() );

	set_a_form_range( option[ OptionKeys::recces::a_form_range ]() );
	set_base_pair_rmsd_cutoff( option[ OptionKeys::recces::base_pair::rmsd_cutoff ]() );
	set_base_pair_rotation_mag( option[ OptionKeys::recces::base_pair::rotation_mag ]() );
	set_base_pair_translation_mag( option[ OptionKeys::recces::base_pair::translation_mag ]() );
	set_sample_jump( option[ OptionKeys::recces::base_pair::sample_jump ]() );

	set_dump_pdb( option[ OptionKeys::recces::dump_pdb ]() );
	set_save_scores( option[ OptionKeys::recces::save_score_terms ]() );
	set_block_stack( option[ OptionKeys::recces::base_pair::block_stack ]() );

	set_dump_freq( option[ OptionKeys::recces::dump_freq ]() );
	set_dump_silent( option[ OptionKeys::recces::dump_silent ]() );
	set_out_torsions( option[ OptionKeys::recces::out_torsions ]() );
	set_setup_base_pair_constraints( option[ OptionKeys::recces::thermal_sampling::setup_base_pair_constraints ]() );

	set_seq1( option[ OptionKeys::recces::seq1 ]() );
	set_seq2( option[ OptionKeys::recces::seq2 ]() );
	if ( seq1().size() > 0 ) {
		legacy_turner_mode_ = true;
	}

	if ( option[ in::file::s ].user() &&
			option[ in::file::s ]().size() > 0 ) {
		set_infile( option[ in::file::s ]()[ 1 ] );
		runtime_assert( !legacy_turner_mode() );
	}
	set_out_prefix( option[ OptionKeys::recces::out_prefix ]() );

	set_temperatures( option[ OptionKeys::recces::temps ]() );
	set_st_weights( option[ OptionKeys::recces::st_weights ]() );

	set_xyz_file( std::string( option[ out::path::path]() ) + "/xyz.txt" );
	set_silent_file( std::string( option[ out::path::path]() ) + option[ out::file::silent ]() );

	if ( option[ OptionKeys::recces::histogram_max ].user() ) set_histogram_max( option[ OptionKeys::recces::histogram_max ]() );
	if ( option[ OptionKeys::recces::histogram_min ].user() ) set_histogram_min( option[ OptionKeys::recces::histogram_min ]() );
	if ( option[ OptionKeys::recces::histogram_spacing ].user() ) set_histogram_spacing( option[ OptionKeys::recces::histogram_spacing ]() );

	// thermal sampling stuff
	set_sample_residues( option[ OptionKeys::recces::thermal_sampling::sample_residues ]() );
	set_free_residues( option[ OptionKeys::recces::thermal_sampling::free_residues ]() );
	set_loop_residues( option[ OptionKeys::recces::thermal_sampling::loop_residues ]() );
	set_angle_range_bb( option[ OptionKeys::recces::thermal_sampling::angle_range_bb ]() );
	set_angle_range_free_bb( option[ OptionKeys::recces::thermal_sampling::angle_range_free_bb ]() );
	set_angle_range_chi( option[ OptionKeys::recces::thermal_sampling::angle_range_chi ]() );
	set_angle_range_free_chi( option[ OptionKeys::recces::thermal_sampling::angle_range_free_chi ]() );
	set_chi_stdev( option[ OptionKeys::recces::thermal_sampling::chi_stdev ]() );
	set_bb_stdev( option[ OptionKeys::recces::thermal_sampling::bb_stdev ]() );
	set_standard_bb_stdev( option[ OptionKeys::recces::thermal_sampling::standard_bb_stdev ]() );
	if ( option[ OptionKeys::recces::thermal_sampling::sample_residues ].user() ) {
		set_thermal_sampler_mode( true );
	}

	// quick hack. probably should change namespace
	rna_secstruct_ = core::pose::rna::RNA_SecStruct( option[ OptionKeys::rna::farna::secstruct ](),
																									 option[ OptionKeys::rna::farna::secstruct_file ]() );
}

} //options
} //recces
} //protocols
