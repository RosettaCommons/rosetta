// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Phil Bradley


// Rosetta headers
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/types.hh>

// ObjexxFCL headers


// Numeric headers


// Utility headers

// C++ headers
#include <map>
#ifdef WIN32
#include <string>
#endif
#include <iostream>

#include <utility/vector1.hh>


namespace core {
namespace scoring {

//////////////////////////////////////////////////////////////////////////////
/// @brief give a ScoreType string name and return its enum type
ScoreType
score_type_from_name( std::string const & name )
{
	return ScoreTypeManager::score_type_from_name( name );
}

//////////////////////////////////////////////////////////////////////////////
/// @brief input operator for ScoreType enum type
///
/// @details read in a string name from a file or std::cin and directly convert
/// it to an ScoreType enum type, for example, std::cin >> ScoreType. This will first check
/// if the lookup map has been set up already. If the string name cannot be
/// converted properly, it will flag the input stream as failure
/// (e.g., istream.fail() is true) and set ScoreType enum type to total_score.
std::istream &
operator >>(
	std::istream & is,
	ScoreType & score_type
)
{
	std::string name;
	is >> name;
	if ( ScoreTypeManager::is_score_type( name ) ) {
		//std::cout << "score_typeextract succeeded " << name << std::endl;
		score_type = ScoreTypeManager::score_type_from_name( name );
	} else {
		std::cout << "score_typeextract failed: " << name << std::endl;
		score_type = total_score;
		is.setstate( std::ios_base::failbit );
	}
	return is;
}

//////////////////////////////////////////////////////////////////////////////
/// @brief output operator for ScoreType enum type
///
/// @details example usage: std::cout << score_type_gly << std::endl;
std::ostream &
operator <<(
	std::ostream & os,
	ScoreType const & score_type
)
{
	os << ScoreTypeManager::name_from_score_type( score_type );
	return os;
}

//////////////////////////////////////////////////////////////////////////////
/// @brief output operator for ScoreTypes list type
std::ostream &
operator <<(
	std::ostream & os,
	ScoreTypes const & score_types
)
{
	for ( core::Size ii(1); ii <= score_types.size(); ++ii ) {
		if ( ii == 1 ) {
			os << "( ";
		} else {
			os << ", ";
		}
		os << score_types[ii];
	}
	os << " )";
	return os;
}

std::string
name_from_score_type( ScoreType  score_type ) {
	return ScoreTypeManager::name_from_score_type( score_type );
}


} // chemical
} // core
