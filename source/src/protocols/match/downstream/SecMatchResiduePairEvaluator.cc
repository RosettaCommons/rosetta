// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/downstream/SecMatchResiduePairEvaluator.cc
/// @brief
/// @author Florian Richter, floric@u.washington.edu, june 09

// Unit headers
#include <protocols/match/downstream/SecMatchResiduePairEvaluator.hh>

// Package headers

// Project headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers

#include <utility/vector1_bool.hh>


namespace protocols {
namespace match {
namespace downstream {

SecMatchResiduePairEvaluator::SecMatchResiduePairEvaluator() = default;


SecMatchResiduePairEvaluator::~SecMatchResiduePairEvaluator() = default;

bool
SecMatchResiduePairEvaluator::require_candidate_residue_atoms_to_lie_near_target_atom( Size /*target_atom_id*/ ) const
{
	return false;
}

utility::vector1< SecMatchResiduePairEvaluator::Size >
SecMatchResiduePairEvaluator::candidate_res_atoms_reqd_near_target_atom(
	Size /*target_atom_id*/
) const {
	// Base class noop.
	utility::vector1< Size > empty;
	return empty;
}

SecMatchResiduePairEvaluator::Real
SecMatchResiduePairEvaluator::max_separation_dist_to_target_atom( Size /*target_atom_id*/ ) const
{
	return -1.0;
}


}
}
}
