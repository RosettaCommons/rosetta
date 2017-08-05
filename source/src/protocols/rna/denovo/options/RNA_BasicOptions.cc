// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/denovo/options/RNA_BasicOptions.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/rna/denovo/options/RNA_BasicOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/OptionKeyList.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.rna.denovo.options.RNA_BasicOptions" );

using namespace basic::options;
using namespace basic::options::OptionKeys;

namespace protocols {
namespace rna {
namespace denovo {
namespace options {

//Constructor
RNA_BasicOptions::RNA_BasicOptions():
	dump_pdb_( false ),
	move_first_rigid_body_( false ),
	verbose_( true )
{}

//Destructor
RNA_BasicOptions::~RNA_BasicOptions()
{}


/// @brief copy constructor
RNA_BasicOptions::RNA_BasicOptions( RNA_BasicOptions const & src ) :
	ResourceOptions( src )
{
	*this = src;
}

/// @brief clone the options
RNA_BasicOptionsOP
RNA_BasicOptions::clone() const
{
	return RNA_BasicOptionsOP( new RNA_BasicOptions( *this ) );
}

///////////////////////////////////////////////////////////////////
void
RNA_BasicOptions::initialize_from_command_line() {
	initialize_from_options( basic::options::option );
}

void
RNA_BasicOptions::initialize_from_options( utility::options::OptionCollection const & opts ) {
	set_dump_pdb( opts[ basic::options::OptionKeys::rna::denovo::out::dump ] ) ;
	set_move_first_rigid_body( opts[ basic::options::OptionKeys::rna::denovo::move_first_rigid_body ] );
}
void
RNA_BasicOptions::list_options_read( utility::options::OptionKeyList & opts ) {
	using namespace basic::options;
	opts + OptionKeys::rna::denovo::out::dump
		+ OptionKeys::rna::denovo::move_first_rigid_body;
}


} //options
} //denovo
} //rna
} //protocols
