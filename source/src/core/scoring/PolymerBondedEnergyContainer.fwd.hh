// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/PolymerBondedEnergyContainer.hh
/// @brief  A container interface long range energies for polymer-bonded residue interactions only.
/// @author Frank DiMaio

#ifndef INCLUDED_core_scoring_PolymerBondedEnergyContainer_fwd_hh
#define INCLUDED_core_scoring_PolymerBondedEnergyContainer_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


namespace core {
namespace scoring {

class PolymerBondedNeighborIterator;
typedef utility::pointer::shared_ptr< PolymerBondedNeighborIterator > PolymerBondedNeighborIteratorOP;
typedef utility::pointer::shared_ptr< const PolymerBondedNeighborIterator > PolymerBondedNeighborIteratorCOP;


class PolymerBondedEnergyContainer;
typedef utility::pointer::shared_ptr< PolymerBondedEnergyContainer > PolymerBondedEnergyContainerOP;
typedef utility::pointer::shared_ptr< const PolymerBondedEnergyContainer > PolymerBondedEnergyContainerCOP;

} // namespace scoring
} // namespace core

#endif
