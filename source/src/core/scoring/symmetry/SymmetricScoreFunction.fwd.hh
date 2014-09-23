// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file   core/scoring/symmetry/SymmetricScoreFunction.fwd.hh
/// @brief  core::scoring::symmetry::SymmetricScoreFunction forward declarations
/// @author Ingemar Andre


#ifndef INCLUDED_core_scoring_symmetry_SymmetricScoreFunction_fwd_hh
#define INCLUDED_core_scoring_symmetry_SymmetricScoreFunction_fwd_hh

#include <utility/pointer/owning_ptr.hh>


namespace core {
namespace scoring {
namespace symmetry {


// Forward
class SymmetricScoreFunction;

typedef utility::pointer::shared_ptr< SymmetricScoreFunction > SymmetricScoreFunctionOP;
typedef utility::pointer::shared_ptr< SymmetricScoreFunction const > SymmetricScoreFunctionCOP;

} // symmetry
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_symmetry_SymmetricScoreFunction_FWD_HH
