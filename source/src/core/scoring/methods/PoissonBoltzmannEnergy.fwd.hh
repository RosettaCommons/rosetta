// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/PoissonBoltzmannEnergy.fwd.hh
/// @brief  Ramachandran energy method class forward declaration
/// @author Phil Bradley
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_scoring_methods_PoissonBoltzmannEnergy_FWD_HH
#define INCLUDED_core_scoring_methods_PoissonBoltzmannEnergy_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace scoring {
namespace methods {

class PoissonBoltzmannEnergy;

typedef utility::pointer::shared_ptr< PoissonBoltzmannEnergy > PoissonBoltzmannEnergyOP;
typedef utility::pointer::shared_ptr< PoissonBoltzmannEnergy const > PoissonBoltzmannEnergyCOP;
typedef utility::pointer::weak_ptr< PoissonBoltzmannEnergy > PoissonBoltzmannEnergyAP;
typedef utility::pointer::weak_ptr< PoissonBoltzmannEnergy const > PoissonBoltzmannEnergyCAP;

class PBLifetimeCache;

typedef utility::pointer::shared_ptr< PBLifetimeCache > PBLifetimeCacheOP;
typedef utility::pointer::shared_ptr< PBLifetimeCache const > PBLifetimeCacheCOP;
typedef utility::pointer::weak_ptr< PBLifetimeCache > PBLifetimeCacheAP;
typedef utility::pointer::weak_ptr< PBLifetimeCache const > PBLifetimeCacheCAP;
} // methods
} // scoring
} // core


#endif
