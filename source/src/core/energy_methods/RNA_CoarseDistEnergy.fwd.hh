// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/energy_methods/RNA_CoarseDistEnergy.fwd.hh
/// @brief Score two-body energies in coarse RNA poses between P, S, and CEN using a statistical potential.
/// @author Andy Watkins (andy.watkins2@gmail.com)

#ifndef INCLUDED_core_energy_methods_RNA_CoarseDistEnergy_fwd_hh
#define INCLUDED_core_energy_methods_RNA_CoarseDistEnergy_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace energy_methods {

class RNA_CoarseDistEnergy;

using RNA_CoarseDistEnergyOP = utility::pointer::shared_ptr< RNA_CoarseDistEnergy >;
using RNA_CoarseDistEnergyCOP = utility::pointer::shared_ptr< RNA_CoarseDistEnergy const >;

} //energy_methods
} //core

#endif //INCLUDED_core_energy_methods_RNA_CoarseDistEnergy_fwd_hh
