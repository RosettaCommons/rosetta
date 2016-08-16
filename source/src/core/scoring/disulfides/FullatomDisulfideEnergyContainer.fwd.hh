// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/disulfides/ConstraintsEnergyContainer.hh
/// @brief  Constraints Energy Container class declaration
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_scoring_disulfides_FullatomDisulfideEnergyContainer_fwd_hh
#define INCLUDED_core_scoring_disulfides_FullatomDisulfideEnergyContainer_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace scoring {
namespace disulfides {

class DisulfResNeighbIterator;

class DisulfResNeighbConstIterator;

class FullatomDisulfideEnergyContainer;

typedef utility::pointer::shared_ptr< FullatomDisulfideEnergyContainer > FullatomDisulfideEnergyContainerOP;
typedef utility::pointer::shared_ptr< FullatomDisulfideEnergyContainer const > FullatomDisulfideEnergyContainerCOP;
typedef utility::pointer::weak_ptr< FullatomDisulfideEnergyContainer > FullatomDisulfideEnergyContainerAP;
typedef utility::pointer::weak_ptr< FullatomDisulfideEnergyContainer const > FullatomDisulfideEnergyContainerCAP;

}
}
}

#endif
