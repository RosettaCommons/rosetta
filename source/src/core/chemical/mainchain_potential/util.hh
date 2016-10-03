// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/mainchain_potential/util.hh
/// @brief  Headers for utility functions for mainchain torsional potentials.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_core_chemical_mainchain_potential_util_hh
#define INCLUDED_core_chemical_mainchain_potential_util_hh

// Unit Headers
#include <core/chemical/mainchain_potential/MainchainScoreTable.fwd.hh>

// Project Headers
//#include <core/types.hh>
//#include <core/pose/Pose.fwd.hh>

// Utility Headers
//#include <utility/pointer/ReferenceCount.hh>
//#include <utility/io/izstream.hh>
#include <utility/vector1.hh>

// Numeric Headers
//#include <numeric/MathNTensor.hh>
//#include <numeric/interpolation/spline/BicubicSpline.hh>

// C++ Headers
#include <map>
#include <sstream>

namespace core {
namespace chemical {
namespace mainchain_potential {

/// @brief Read a Shapovalov/Ramachandran-style mainchain torsion file, parse it for mainchain
/// potentials corresponding to a vector of ResidueType names, and return a map of (name->MainchainScoreTableCOP).
/// @param[in] filename The name of the file to read.
/// @param[in] res_type_names A vector of ResidueType names that this function will seek data for.
/// @param[in] use_polycubic_interpolation Should the MainchainScoreTables use polycubic interpolation?
/// @param[out] mainchain_score_table_map The output map of (name->MainchainScoreTableCOP).
void read_rama_map_file_shapovalov(
	std::string const &filename,
	utility::vector1< std::string > const &res_type_names,
	bool const use_polycubic_interpolation,
	std::map< std::string, MainchainScoreTableCOP > &mainchain_score_table_map
);

/// @brief Read a Shapovalov/Ramachandran-style mainchain torsion file, parse it for all the mainchain
/// potentials that it contains, and return a vector of pairs of (map name, MainchainScoreTableOP).
/// @param[in] filename The name of the file to read.
/// @param[in] use_polycubic_interpolation Should the MainchainScoreTables use polycubic interpolation?
/// @param[out] newtables The output vector of pairs of (name, MainchainScoreTableOP).
void read_rama_map_file_shapovalov(
	std::string const &filename,
	bool const use_polycubic_interpolation,
	utility::vector1< std::pair < std::string, MainchainScoreTableOP> > &newtables
);

} //mainchain_potential
} //chemical
} //core

#endif //INCLUDED_core_chemical_mainchain_potential_util_hh
