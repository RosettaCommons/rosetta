// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/UnfoldedStatePotential.fwd.hh
/// @brief  Unfolded state energies based on energies of residues in fragments, forward declaration
/// @author Ron Jacak (ronj@email.unc.edu)

#ifndef INCLUDED_core_scoring_UnfoldedStatePotential_fwd_hh
#define INCLUDED_core_scoring_UnfoldedStatePotential_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {

class UnfoldedStatePotential;

typedef utility::pointer::shared_ptr< UnfoldedStatePotential > UnfoldedStatePotentialOP;

} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_UnfoldedStatePotential_FWD_HH
