// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/symmetry/SymmetricScoreFunction.hh
/// @brief  Symmetric Score function class
/// @author Ingemar Andre


#ifndef INCLUDED_core_scoring_symmetry_SymmetricScoreFunction_hh
#define INCLUDED_core_scoring_symmetry_SymmetricScoreFunction_hh

// Unit headers
#include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>

// Package headers
#include <core/scoring/ScoreFunction.hh>

// Project headers
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>

#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace symmetry {

/// @brief This class is provided as a backward compatibility shim
/// At this point it should be identical to the base ScoreFunction.
class SymmetricScoreFunction : public ScoreFunction
{
	using ScoreFunction::ScoreFunction;
};

} // symmetry
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_SymmetricScoreFunction_HH
