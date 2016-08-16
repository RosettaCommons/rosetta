// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ScoreTypeManager.hh
/// @brief  Score type enumeration
/// @author Phil Bradley


#ifndef INCLUDED_core_scoring_ScoreTypeManager_hh
#define INCLUDED_core_scoring_ScoreTypeManager_hh

// Project Headers
#include <core/scoring/ScoreType.hh>

// Utility Headers
#include <utility/vector1.hh>

// STL Headers
#include <map>


namespace core {
namespace scoring {

class ScoreTypeManager {
public:
	static
	ScoreType
	score_type_from_name( std::string const & name );

	static
	std::string
	name_from_score_type( ScoreType score_type );

	static
	bool
	is_score_type( std::string const & name );

private:
	static void setup_score_type_names();

private:
	static bool initialized_;

	// lookup map from string name to enum type
	static std::map< std::string, ScoreType > name2score_type_;
	static utility::vector1< std::string >    score_type2name_;

}; // ScoreTypeManager

} // scoring
} // core

#endif
