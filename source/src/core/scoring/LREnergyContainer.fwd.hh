// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/LongRangeEnergyContainer.fwd.hh
/// @brief  A container for storing and scoring long range energies forward declaration
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_scoring_LREnergyContainer_fwd_hh
#define INCLUDED_core_scoring_LREnergyContainer_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {

class ResidueNeighborIterator;
typedef utility::pointer::shared_ptr< ResidueNeighborIterator > ResidueNeighborIteratorOP;
typedef utility::pointer::shared_ptr< ResidueNeighborIterator const > ResidueNeighborIteratorCOP;

class ResidueNeighborConstIterator;
typedef utility::pointer::shared_ptr< ResidueNeighborConstIterator > ResidueNeighborConstIteratorOP;
typedef utility::pointer::shared_ptr< ResidueNeighborConstIterator const > ResidueNeighborConstIteratorCOP;


class LREnergyContainer;
typedef utility::pointer::shared_ptr< LREnergyContainer > LREnergyContainerOP;
typedef utility::pointer::shared_ptr< LREnergyContainer const > LREnergyContainerCOP;

}
}

#endif
