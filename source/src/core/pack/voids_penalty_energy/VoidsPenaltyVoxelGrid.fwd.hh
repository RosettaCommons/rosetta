// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/voids_penalty_energy/VoidsPenaltyVoxelGrid.fwd.hh
/// @brief A 3D boolean array used for identifying core voxels in the VoidsPenaltyEnergy.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_core_pack_voids_penalty_energy_VoidsPenaltyVoxelGrid_fwd_hh
#define INCLUDED_core_pack_voids_penalty_energy_VoidsPenaltyVoxelGrid_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace pack {
namespace voids_penalty_energy {

class VoidsPenaltyVoxelGrid;

typedef utility::pointer::shared_ptr< VoidsPenaltyVoxelGrid > VoidsPenaltyVoxelGridOP;
typedef utility::pointer::shared_ptr< VoidsPenaltyVoxelGrid const > VoidsPenaltyVoxelGridCOP;

} //core
} //pack
} //voids_penalty_energy

#endif //INCLUDED_core_pack_voids_penalty_energy_VoidsPenaltyVoxelGrid_fwd_hh
