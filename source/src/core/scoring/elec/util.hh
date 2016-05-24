// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/util.hh
/// @brief  Utility functions for the distance-dependent dielectric potential.
/// @author Vikram K. Mulligan (vmullig@uw.edu).


#ifndef INCLUDED_core_scoring_elec_util_hh
#define INCLUDED_core_scoring_elec_util_hh

// Unit Headers
#include <core/scoring/elec/FA_ElecEnergy.fwd.hh>
#include <core/scoring/elec/CPRepMapType.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ headers


namespace core {
namespace scoring {
namespace elec {

/// @brief Read the CP tables from the database and return an owning pointer to the
/// new object created in memory.
/// @details Called by the ScoringManager to allow these data to be read in once and
/// only once.
CPRepMapTypeOP read_cp_tables_from_db( std::string const &filename );

} // namespace elec
} // namespace scoring
} // namespace core

#endif
