// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/import_pose/options/RNA_DeNovoProtocolOptions.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/import_pose/options/RNA_DeNovoProtocolOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/OptionKeyList.hh>


static basic::Tracer TR( "core.import_pose.options.RNA_DeNovoProtocolOptions" );

using namespace basic::options;
using namespace basic::options::OptionKeys;

namespace core {
namespace import_pose {
namespace options {

//Constructor
RNA_DeNovoProtocolOptions::RNA_DeNovoProtocolOptions():
	nstruct_( 1 ),
	lores_scorefxn_( "rna/denovo/rna_lores.wts" ),
	align_output_( true ),
	overwrite_( false ),
	binary_rna_output_( false ),
	use_legacy_setup_( false ),
	cst_gap_( false ),
	dump_stems_( false )
{}

//Destructor
RNA_DeNovoProtocolOptions::~RNA_DeNovoProtocolOptions() = default;

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
	RNA_FragmentMonteCarloOptions::initialize_from_options( basic::options::option );
	initialize_from_options( basic::options::option );
}

void
RNA_DeNovoProtocolOptions::initialize_from_options( utility::options::OptionCollection const & opts ) {
	RNA_FragmentMonteCarloOptions::initialize_from_options( opts );

	nstruct_ = opts[ out::nstruct ]();
	if ( opts[ out::file::silent ].user() ) {
		silent_file_ = opts[ out::file::silent ]();
	} else {
		std::string tag;
		if ( opts[ basic::options::OptionKeys::rna::denovo::tag ].user() ) {
			tag = opts[ basic::options::OptionKeys::rna::denovo::tag ]();
		} else {
			tag = utility::file_basename( utility::file::cwd() );
			TR << TR.Green << "Setting silent file name based on directory: " << tag << ".out" << std::endl;
		}
		silent_file_ = tag  + ".out";
	}

	set_overwrite(  opts[ out::overwrite ] );

	if ( option[ basic::options::OptionKeys::edensity::mapfile ].user() ) {
		// should we have an option that allows this to be
		// turned on again if we have density (?)
		set_align_output( false );
	}

	// note that althrough the following variables are held in the base class RNA_FragmentMonteCarloOptions, they are not initialized from command-line there.
	// they really should only be set up for runs using the rna_denovo exectuable -- so they are set up here.
	if ( opts[ in::file::s ].user() ) set_chunk_pdb_files( opts[ in::file::s ]() );
	if ( opts[ basic::options::OptionKeys::rna::denovo::initial_structures ].user() ) {
		set_chunk_initialization_pdb_files( opts[ basic::options::OptionKeys::rna::denovo::initial_structures ]() );
	}
	if ( opts[ in::file::silent ].user() )  set_chunk_silent_files( opts[ in::file::silent ]() );
	if ( opts[ in::file::input_res ].user() )  set_input_res( opts[ in::file::input_res ]() ) ;

	if ( opts[ basic::options::OptionKeys::rna::denovo::lores_scorefxn ].user() ) set_lores_scorefxn( opts[ basic::options::OptionKeys::rna::denovo::lores_scorefxn ] );

	if ( opts[ in::file::silent_struct_type ]() == "binary_rna"  ||
			opts[ basic::options::OptionKeys::rna::denovo::out::binary_output ]() ||
			close_loops() ||
			vary_bond_geometry() ) set_binary_rna_output( true );

	use_legacy_setup_ = opts[ basic::options::OptionKeys::rna::denovo::use_legacy_setup ]();

	cst_gap_ = opts[ basic::options::OptionKeys::rna::denovo::cst_gap ]();

	dump_stems_ = option[ basic::options::OptionKeys::rna::denovo::dump_stems ]();
}

void
RNA_DeNovoProtocolOptions::list_options_read( utility::options::OptionKeyList & opts ) {
	using namespace basic::options;
	RNA_FragmentMonteCarloOptions::list_options_read( opts );
	opts + OptionKeys::out::nstruct
		+ OptionKeys::out::file::silent
		+ OptionKeys::rna::denovo::tag
		+ OptionKeys::rna::denovo::out::output_lores_silent_file
		+ OptionKeys::out::overwrite
		+ OptionKeys::in::file::s
		+ OptionKeys::in::file::silent
		+ OptionKeys::in::file::input_res
		+ OptionKeys::rna::denovo::lores_scorefxn
		+ OptionKeys::in::file::silent_struct_type
		+ OptionKeys::rna::denovo::out::binary_output
		+ OptionKeys::rna::denovo::use_legacy_setup
		+ OptionKeys::rna::denovo::cst_gap
		+ OptionKeys::rna::denovo::dump_stems
		+ OptionKeys::rna::denovo::initial_structures;
}

} //denovo
} //rna
} //protocols
