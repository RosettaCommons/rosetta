// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/carbohydrates/database_io.hh
/// @brief   Database input/output function declarations for carbohydrate-specific scoring data.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_scoring_carbohydrates_database_io_HH
#define INCLUDED_core_scoring_carbohydrates_database_io_HH

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ headers
#include <map>
#include <string>


namespace core {
namespace scoring {
namespace carbohydrates {

/// @brief  Return a table of Gaussian parameters read from a database file.
std::map< char, utility::vector1< Real > > read_Gaussian_parameters_from_database_file( std::string const & filename );

}  // namespace carbohydrates
}  // namespace scoring
}  // namespace core

#endif  // INCLUDED_core_scoring_carbohydrates_database_io_HH
