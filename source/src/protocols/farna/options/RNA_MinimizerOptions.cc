// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/farna/options/RNA_MinimizerOptions.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/farna/options/RNA_MinimizerOptions.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.farna.options.RNA_MinimizerOptions" );

using namespace basic::options;
using namespace basic::options::OptionKeys;

namespace protocols {
namespace farna {
namespace options {

	//Constructor
	RNA_MinimizerOptions::RNA_MinimizerOptions():
		minimize_rounds_( 2 ),
		deriv_check_( false ),
		skip_o2prime_trials_( false ),
		vary_bond_geometry_( false ),
    minimizer_use_coordinate_constraints_( false ),
    minimize_bps_( false )
	{}

	//Destructor
	RNA_MinimizerOptions::~RNA_MinimizerOptions()
	{}


	/// @brief copy constructor
	RNA_MinimizerOptions::RNA_MinimizerOptions( RNA_MinimizerOptions const & src ) :
		ResourceOptions( src ),
		RNA_BasicOptions( src )
	{
		*this = src;
	}

	/// @brief clone the options
	RNA_MinimizerOptionsOP
	RNA_MinimizerOptions::clone() const
	{
		return RNA_MinimizerOptionsOP( new RNA_MinimizerOptions( *this ) );
	}

	///////////////////////////////////////////////////////////////////
	void
	RNA_MinimizerOptions::initialize_from_command_line() {

		RNA_BasicOptions::initialize_from_command_line();

		set_minimize_rounds( option[ OptionKeys::rna::farna::minimize::minimize_rounds ]() );
		set_vary_bond_geometry( option[ OptionKeys::rna::vary_geometry ]() );
		set_deriv_check( option[ OptionKeys::rna::farna::minimize::deriv_check ]() );
		set_skip_o2prime_trials( option[ OptionKeys::rna::farna::minimize::skip_o2prime_trials]() );

		set_extra_minimize_res( option[ OptionKeys::rna::farna::minimize::extra_minimize_res ]() ) ;
		set_extra_minimize_chi_res( option[ OptionKeys::rna::farna::minimize::extra_minimize_chi_res ]() ) ;

		set_minimizer_use_coordinate_constraints( option[ OptionKeys::rna::farna::minimize::minimizer_use_coordinate_constraints ]() );
		if ( option[ OptionKeys::rna::farna::minimize::skip_coord_constraints]() ) set_minimizer_use_coordinate_constraints( false );
		set_minimize_bps( option[ OptionKeys::rna::farna::minimize::minimize_bps ]() ) ;

	}

} //options
} //farna
} //protocols
