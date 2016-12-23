// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/farna/options/RNA_DeNovoProtocolOptions.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/farna/options/RNA_DeNovoProtocolOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/file/file_sys_util.hh>


static basic::Tracer TR( "protocols.farna.options.RNA_DeNovoProtocolOptions" );

using namespace basic::options;
using namespace basic::options::OptionKeys;

namespace protocols {
namespace farna {
namespace options {

//Constructor
RNA_DeNovoProtocolOptions::RNA_DeNovoProtocolOptions():
	nstruct_( 1 ),
	lores_scorefxn_( "farna/rna_lores.wts" ),
	output_lores_silent_file_( false ),
	output_filters_( false ),
	binary_rna_output_( false ),
	save_times_( false ),
	use_legacy_setup_( false ),
	cst_gap_( false )
{}

//Destructor
RNA_DeNovoProtocolOptions::~RNA_DeNovoProtocolOptions()
{}

/// @brief copy constructor
RNA_DeNovoProtocolOptions::RNA_DeNovoProtocolOptions( RNA_DeNovoProtocolOptions const & src ) :
	ResourceOptions( src ),
	RNA_BasicOptions( src ),
	RNA_MinimizerOptions( src ),
	RNA_FragmentMonteCarloOptions( src )
{
	*this = src;
}

/// @brief clone the options
RNA_DeNovoProtocolOptionsOP
RNA_DeNovoProtocolOptions::clone() const
{
	return RNA_DeNovoProtocolOptionsOP( new RNA_DeNovoProtocolOptions( *this ) );
}

///////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocolOptions::initialize_from_command_line() {

	RNA_FragmentMonteCarloOptions::initialize_from_command_line();

	nstruct_ = option[ out::nstruct ]();
	if ( option[ out::file::silent ].user() ) {
		silent_file_ = option[ out::file::silent ]();
	} else {
		std::string tag;
		if ( option[ rna::farna::tag ].user() ) {
			tag = option[ rna::farna::tag ]();
		} else {
			tag = utility::file_basename( utility::file::cwd() );
			TR << TR.Green << "Setting silent file name based on directory: " << tag << ".out" << std::endl;
		}
		silent_file_ = tag  + ".out";
	}

	set_output_lores_silent_file( option[ rna::farna::output_lores_silent_file ] );
	set_output_filters(  option[ rna::farna::output_filters ] );

	// note that althrough the following variables are held in the base class RNA_FragmentMonteCarloOptions, they are not initialized from command-line there.
	// they really should only be set up for runs using the rna_denovo exectuable -- so they are set up here.
	if ( option[ in::file::s ].user() ) set_chunk_pdb_files( option[ in::file::s ]() );
	if ( option[ in::file::silent ].user() )  set_chunk_silent_files( option[ in::file::silent ]() );
	if ( option[ in::file::input_res ].user() )  set_input_res( option[ in::file::input_res ]() ) ;

	if ( option[ rna::farna::lores_scorefxn ].user() ) set_lores_scorefxn( option[ rna::farna::lores_scorefxn ] );

	if ( option[ in::file::silent_struct_type ]() == "binary_rna"  ||
			option[ rna::farna::binary_output ]() ||
			close_loops() ||
			vary_bond_geometry() ) set_binary_rna_output( true );

	save_times_ = option[ OptionKeys::out::save_times ]();

	use_legacy_setup_ = option[ rna::farna::use_legacy_setup ]();

	cst_gap_ = option[ rna::farna::cst_gap ]();
}

} //options
} //farna
} //protocols
