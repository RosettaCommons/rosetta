// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/PeptideBondedEnergyContainer.hh
/// @brief  A container interface long range energies for n->n+1 interactions only
/// @author Frank DiMaio

#ifndef INCLUDED_core_scoring_PeptideBondedEnergyContainer_fwd_hh
#define INCLUDED_core_scoring_PeptideBondedEnergyContainer_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


namespace core {
namespace scoring {

class PeptideBondedNeighborIterator;
typedef utility::pointer::shared_ptr< PeptideBondedNeighborIterator > PeptideBondedNeighborIteratorOP;
typedef utility::pointer::shared_ptr< const PeptideBondedNeighborIterator > PeptideBondedNeighborIteratorCOP;


class PeptideBondedEnergyContainer;
typedef utility::pointer::shared_ptr< PeptideBondedEnergyContainer > PeptideBondedEnergyContainerOP;
typedef utility::pointer::shared_ptr< const PeptideBondedEnergyContainer > PeptideBondedEnergyContainerCOP;

} // namespace scoring
} // namespace core

#endif
