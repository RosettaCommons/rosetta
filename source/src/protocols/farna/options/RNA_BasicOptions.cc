// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/farna/options/RNA_BasicOptions.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/farna/options/RNA_BasicOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.farna.options.RNA_BasicOptions" );

using namespace basic::options;
using namespace basic::options::OptionKeys;

namespace protocols {
namespace farna {
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
		set_dump_pdb( option[ rna::farna::dump ] ) ;
		set_move_first_rigid_body(  option[ rna::farna::move_first_rigid_body ] );

	}

} //options
} //farna
} //protocols
