// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/etable/count_pair/CountPairCrossover3.hh
/// @brief  Count pair for residue pairs connected with one bond, where the
/// crossover from excluding to counting atom pair interactions is at 3 bonds.
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


// Unit Headers
#include <core/scoring/etable/count_pair/CountPairCrossover3.hh>

// Package Headers
#include <core/scoring/etable/count_pair/CountPairFunction.hh>

// Project Headers
#include <core/types.hh>

namespace core {
namespace scoring {
namespace etable {
namespace count_pair {

// @brief virtual destrutor
CountPairCrossover3::~CountPairCrossover3() = default;


} // namespace count_pair
} // namespace etable
} // namespace scoring
} // namespace core

