// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @brief Forward header file for class ExactOccludedHbondSolEnergy
///
/// @author Andrea Bazzoli (bazzoli@ku.edu)


#ifndef INCLUDED_ExactOccludedHbondSolEnergyCOP_fwd_hh
#define INCLUDED_ExactOccludedHbondSolEnergyCOP_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace geometric_solvation {


class ExactOccludedHbondSolEnergy;

typedef  utility::pointer::shared_ptr< ExactOccludedHbondSolEnergy > ExactOccludedHbondSolEnergyOP;
typedef  utility::pointer::shared_ptr< ExactOccludedHbondSolEnergy const> ExactOccludedHbondSolEnergyCOP;

} // geometric_solvation
} // scoring
} // core


#endif
