// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/replica_docking/TempWeightedMetropolisHastingsMover.fwd.hh
/// @brief  TempWeightedMetropolisHastingsMover forward declarations header
/// @author  Zhe Zhang

#ifndef INCLUDED_devel_replica_docking_TempWeightedMetropolisHastingsMover_fwd_hh
#define INCLUDED_devel_replica_docking_TempWeightedMetropolisHastingsMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace devel {
namespace replica_docking {

//Forwards and OP typedefs
class TempWeightedMetropolisHastingsMover;
typedef utility::pointer::owning_ptr< TempWeightedMetropolisHastingsMover > TempWeightedMetropolisHastingsMoverOP;
typedef utility::pointer::owning_ptr< TempWeightedMetropolisHastingsMover const > TempWeightedMetropolisHastingsMoverCOP;
typedef utility::pointer::access_ptr< TempWeightedMetropolisHastingsMover > TempWeightedMetropolisHastingsMoverAP;
typedef utility::pointer::access_ptr< TempWeightedMetropolisHastingsMover const > TempWeightedMetropolisHastingsMoverCAP;

} //
} //

#endif //INCLUDED_devel_replica_docking_TempWeightedMetropolisHastingsMover_FWD_HH
