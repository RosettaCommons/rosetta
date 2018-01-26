// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/import_pose/options/RNA_BasicOptions.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/import_pose/options/RNA_BasicOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/OptionKeyList.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "core.import_pose.options.RNA_BasicOptions" );

using namespace basic::options;
using namespace basic::options::OptionKeys;

namespace core {
namespace import_pose {
namespace options {

//Constructor
RNA_BasicOptions::RNA_BasicOptions():
	dump_pdb_( false ),
	move_first_rigid_body_( false ),
	dock_into_density_( false ),
	model_with_density_( false ),
	verbose_( true )
{}

//Destructor
RNA_BasicOptions::~RNA_BasicOptions() = default;


/// @brief copy constructor
RNA_BasicOptions::RNA_BasicOptions( RNA_BasicOptions const & src ) :
	ReferenceCount()
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
	if ( opts[ basic::options::OptionKeys::edensity::mapfile ].user() ) {
		// default false, can only be true if we really have a density map
		set_dock_into_density( option[ basic::options::OptionKeys::rna::denovo::dock_into_density ] );
		set_model_with_density( true );
	}

	set_dump_pdb( opts[ basic::options::OptionKeys::rna::denovo::out::dump ] ) ;
	set_move_first_rigid_body( opts[ basic::options::OptionKeys::rna::denovo::move_first_rigid_body ] );
}

void
RNA_BasicOptions::list_options_read( utility::options::OptionKeyList & opts ) {
	using namespace basic::options;
	opts + OptionKeys::rna::denovo::out::dump
		+ OptionKeys::rna::denovo::move_first_rigid_body
		+ OptionKeys::edensity::mapfile
		+ OptionKeys::rna::denovo::dock_into_density;
}


} //options
} //rna
} //protocols
