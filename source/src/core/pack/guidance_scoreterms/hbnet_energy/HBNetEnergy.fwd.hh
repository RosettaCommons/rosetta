// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/guidance_scoreterms/hbnet_energy/HBNetEnergy.fwd.hh
/// @brief Forward declarations for an EnergyMethod that gives a bonus for hydrogen bond networks, which ramps nonlinearly with the size of the
/// networks.
/// @details This energy method is inherently not pairwise decomposible.  However, it is intended for very rapid calculation,
/// and has been designed to plug into Alex Ford's modifications to the packer that permit it to work with non-pairwise scoring
/// terms.  It has also been modified to permit sub-regions of a pose to be scored.
/// @author Vikram K. Mulligan (vmullig@uw.edu).


#ifndef INCLUDED_core_pack_guidance_scoreterms_hbnet_energy_HBNetEnergy_fwd_hh
#define INCLUDED_core_pack_guidance_scoreterms_hbnet_energy_HBNetEnergy_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace hbnet_energy {

class HBNetEnergy;

typedef utility::pointer::shared_ptr< HBNetEnergy > HBNetEnergyOP;
typedef utility::pointer::shared_ptr< HBNetEnergy const > HBNetEnergyCOP;


} // hbnet_energy
} // guidance_scoreterms
} // pack
} // core


#endif
