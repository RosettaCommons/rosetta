// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/energy_methods/ArgCationPiEnergy.fwd.hh
/// @brief  Cation pi term that specializes in bringing Arginine and rings together
/// @details Currently designed for canonical amino acids but easily extended beyond.
/// @author Brian Coventry (bcov@uw.edu)

#ifndef INCLUDED_core_energy_methods_ArgCationPiEnergy_fwd_hh
#define INCLUDED_core_energy_methods_ArgCationPiEnergy_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace methods {


class ArgCationPiEnergy;

typedef utility::pointer::shared_ptr< ArgCationPiEnergy > ArgCationPiEnergyOP;
typedef utility::pointer::shared_ptr< ArgCationPiEnergy const > ArgCationPiEnergyCOP;


} // methods
} // scoring
} // core

#endif
