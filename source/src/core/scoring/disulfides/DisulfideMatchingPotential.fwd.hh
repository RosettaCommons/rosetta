// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/disulfides/DisulfideMatchingPotential.fwd.hh
/// @author rvernon@u.washington.edu
/// @date   02/09/10


#ifndef INCLUDED_core_scoring_disulfides_DisulfideMatchingPotential_fwd_hh
#define INCLUDED_core_scoring_disulfides_DisulfideMatchingPotential_fwd_hh

//Utility Headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace scoring {
namespace disulfides {

class DisulfideMatchingPotential;
typedef utility::pointer::shared_ptr< DisulfideMatchingPotential > DisulfideMatchingPotentialOP;
typedef utility::pointer::shared_ptr< DisulfideMatchingPotential const > DisulfideMatchingPotentialCOP;
typedef utility::pointer::weak_ptr< DisulfideMatchingPotential > DisulfideMatchingPotentialAP;
typedef utility::pointer::weak_ptr< DisulfideMatchingPotential const > DisulfideMatchingPotentialCAP;


} //disulfides
} //scoring
} //core


#endif //INCLUDED_core_scoring_disulfides_DisulfideMatchingPotential_FWD_HH
