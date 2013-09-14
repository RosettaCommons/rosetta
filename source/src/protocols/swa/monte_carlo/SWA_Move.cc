// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/swa/monte_carlo/SWA_Move.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/swa/monte_carlo/SWA_Move.hh>
#include <map>
#include <iostream>
#include <utility/string_util.hh>

using utility::make_tag_with_dashes;

namespace protocols {
namespace swa {
namespace monte_carlo {

	//Constructor
	SWA_Move::SWA_Move( Chunk const & chunk,
								MovingResidueCase const & moving_residue_case,
								AddOrDeleteChoice const & add_or_delete_choice ):
		utility::pointer::ReferenceCount(),
		chunk_( chunk ),
		moving_residue_case_( moving_residue_case ),
		add_or_delete_choice_( add_or_delete_choice )
	{
	}

	SWA_Move::SWA_Move():
		utility::pointer::ReferenceCount(),
		moving_residue_case_( NO_CASE ),
		add_or_delete_choice_( NO_ADD_OR_DELETE )
	{
	}


	//Destructor
	SWA_Move::~SWA_Move()
	{}

	SWA_Move::SWA_Move( SWA_Move const & src ):
		ReferenceCount( src )
	{
		*this = src;
	}

	SWA_Move &
	SWA_Move::operator=( SWA_Move const & src )
	{
		chunk_ = src.chunk_;
		moving_residue_case_ = src.moving_residue_case_;
		add_or_delete_choice_ = src.add_or_delete_choice_;
		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	std::string
	to_string( MovingResidueCase const & moving_residue_case ){

		static bool init( false );
		static std::map< MovingResidueCase, std::string> moving_residue_case_name;

		if ( !init ){
			moving_residue_case_name[ NO_CASE ] = "NO_CASE";
			moving_residue_case_name[ CHAIN_TERMINUS_5PRIME ] = "CHAIN_TERMINUS_5PRIME";
			moving_residue_case_name[ CHAIN_TERMINUS_3PRIME ] = "CHAIN_TERMINUS_3PRIME";
			moving_residue_case_name[ INTERNAL ] = "INTERNAL";
			moving_residue_case_name[ FLOATING_BASE ] = "FLOATING_BASE";
			init = true;
		}

		return moving_residue_case_name[ moving_residue_case ];
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	std::string
	to_string( AddOrDeleteChoice const & add_or_delete_choice ){

		static bool init( false );
		static std::map< AddOrDeleteChoice, std::string> add_or_delete_choice_name;

		if ( !init ){
			add_or_delete_choice_name[ NO_ADD_OR_DELETE ] = "NO_ADD_OR_DELETE";
			add_or_delete_choice_name[ ADD ]    = "ADD";
			add_or_delete_choice_name[ DELETE ] = "DELETE";
			init = true;
		}

		return add_or_delete_choice_name[ add_or_delete_choice ];
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	std::ostream &
	operator <<( std::ostream & os, SWA_Move const swa_move )
	{
		os << "Choice " << to_string( swa_move.add_or_delete_choice() )	<< "   Case " << to_string( swa_move.moving_residue_case() ) <<
			"   Residues " << make_tag_with_dashes( swa_move.chunk() );
		return os;
	}

} //monte_carlo
} //swa
} //protocols
