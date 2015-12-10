// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/full_model_info/SubMotifInfo.cc
/// @brief  Stores information about submotifs in a pose.
/// @author Caleb Geniesse

// Unit headers
#include <core/pose/full_model_info/SubMotifInfo.hh>

// Package headers

// Project headers

// Utility headers
#include <utility/vector1.hh>

// C++
#include <string>

//////////////////////////////////////////////////////////////////////////////
//
// This object holds information on the 'submotifs' added during stepwise
// modeling. It is used in stepwise monte carlo to determine what moves are
// allowed next for the pose. The SubMotifInfo object is critical for
// satisfying detailed balance when submotif_moves are allowed.
//
//  -- Caleb, 2015
//
// A little more detail:
//
//////////////////////////////////////////////////////////////////////////////

using namespace core;

namespace core {
namespace pose {
namespace full_model_info {


//////////////////////////////////////////////////////////////////////////////
// @brief Constructor
SubMotifInfo::SubMotifInfo() :
	tag_( "" ),
	seed_( false )
{}

//////////////////////////////////////////////////////////////////////////////
// @brief Constructor
SubMotifInfo::SubMotifInfo(
	utility::vector1< Size > const & res_list,
	std::string const & tag,
	bool const & seed /*= false*/
) :
	res_list_( res_list ),
	tag_( tag ),
	seed_( seed )
{}

//////////////////////////////////////////////////////////////////////////////
// @brief Copy Constructor
SubMotifInfo::SubMotifInfo( SubMotifInfo const & src ) :
	res_list_( src.res_list_ ),
	tag_( src.tag_ ),
	seed_( src.seed_ )
{}

//////////////////////////////////////////////////////////////////////////////
// @brief Destructor
SubMotifInfo::~SubMotifInfo()
{}


utility::vector1< Size >
SubMotifInfo::sorted_res_list() const
{
	utility::vector1< Size > sorted_res_list = res_list_;
	std::sort( sorted_res_list.begin(), sorted_res_list.end() );
	return sorted_res_list;
}


/// @brief Equality comparator
bool
operator==( SubMotifInfoOP lhs, SubMotifInfoOP rhs )
{
	if ( !lhs || !rhs ) return false;
	return ( lhs->sorted_res_list() == rhs->sorted_res_list() &&
		lhs->tag() == rhs->tag() &&
		lhs->seed() == rhs->seed() );
}

/// @brief Equality comparator
bool
operator!=( SubMotifInfoOP lhs, SubMotifInfoOP rhs )
{
	return !( lhs == rhs );
}

/// @brief << operator
std::ostream &
operator <<( std::ostream & os, SubMotifInfoOP submotif_info )
{
	os << "SUBMOTIF_INFO";
	os << " RES_LIST";
	for ( Size i = 1; i <= submotif_info->res_list().size(); ++i ) {
		os << " " << submotif_info->sorted_res_list( i );
	}
	os << " SUBMOTIF_TAG" << " " << submotif_info->tag();
	os << " SEED" << " " << submotif_info->seed();
	return os;
}

/// @brief << operator
std::istream &
operator >>( std::istream & is, SubMotifInfoOP submotif_info )
{
	using namespace utility;
	std::string str;

	//
	// SUBMOTIF_INFO RES_LIST 5 3 SUBMOTIF_TAG base_pairs/gg_2grb_RNA.pdb SEED 0
	//

	is >> str;
	runtime_assert( !is.fail() && str == "SUBMOTIF_INFO" );

	is >> str;
	runtime_assert( !is.fail() && str == "RES_LIST" );
	utility::vector1<Size> res_list;
	while ( true ) {
		is >> str;
		if ( str == "SUBMOTIF_TAG" ) {
			break;
		}
		res_list.push_back( atoi(str.c_str()) );
	}
	submotif_info->res_list( res_list );

	runtime_assert( !is.fail() && str == "SUBMOTIF_TAG" );
	is >> str;
	submotif_info->tag( str );

	is >> str;
	runtime_assert( !is.fail() && str == "SEED" );
	submotif_info->seed( atoi(str.c_str()) );

	// decoy tag?
	is >> str;
	return is;
}

} //full_model_info
} //pose
} //core
