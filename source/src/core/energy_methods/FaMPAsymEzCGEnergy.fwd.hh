// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/energy_methods/FaMPAsymEzCG.fwd.hh
///
/// @brief  Fullatom asymetric EZ potential for CG atoms
/// @details Asymetric EZ potential for CG atoms, from Schramm et al 2012 Structure
///
/// @author  Meghan Franklin (meghanwfranklin@gmail.com)

#ifndef INCLUDED_core_energy_methods_FaMPAsymEzCGEnergy_fwd_hh
#define INCLUDED_core_energy_methods_FaMPAsymEzCGEnergy_fwd_hh

// Utility Methods
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace energy_methods {

class FaMPAsymEzCGEnergy;
typedef utility::pointer::shared_ptr< FaMPAsymEzCGEnergy > FaMPAsymEzCGEnergyOP;
typedef utility::pointer::shared_ptr< FaMPAsymEzCGEnergy const > FaMPAsymEzCGEnergyCOP;

} // energy_methods
} // core

#endif // INCLUDED_core_energy_methods_FaMPAsymEzCGEnergy_fwd_hh
