// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/P_AA_Energy.fwd.hh
/// @brief  Probability of observing an amino acid (NOT conditional on phi/psi), energy method forward declaration
/// @author Ron Jacak


#ifndef INCLUDED_core_scoring_methods_P_AA_Energy_fwd_hh
#define INCLUDED_core_scoring_methods_P_AA_Energy_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace methods {

class P_AA_Energy;

typedef utility::pointer::shared_ptr< P_AA_Energy > P_AA_EnergyOP;

} // methods
} // scoring
} // core


#endif // INCLUDED_core_scoring_methods_P_AA_Energy_FWD_HH
