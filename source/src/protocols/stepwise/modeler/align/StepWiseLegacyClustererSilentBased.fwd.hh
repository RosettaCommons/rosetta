// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/Pose.fwd.hh
/// @brief  Various OPS for various classes.
/// @author Rhiju Das

#include <utility/pointer/owning_ptr.hh>

#ifndef INCLUDED_protocols_STEPWISE_Clusterer_FWD_HH
#define INCLUDED_protocols_STEPWISE_Clusterer_FWD_HH

namespace protocols {
namespace stepwise {
namespace modeler {
namespace align {

class StepWiseLegacyClustererSilentBased;
typedef utility::pointer::shared_ptr< StepWiseLegacyClustererSilentBased > StepWiseLegacyClustererSilentBasedOP;
typedef utility::pointer::shared_ptr< StepWiseLegacyClustererSilentBased const > StepWiseLegacyClustererSilentBasedCOP;

} //align
} //modeler
} //stepwise
} //protocols

#endif
