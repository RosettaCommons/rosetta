// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/canonical_sampling/MetropolisHastingsMover.fwd.hh
/// @brief  MetropolisHastingsMover forward declarations header
/// @author

#ifndef INCLUDED_protocols_canonical_sampling_SidechainMetropolisHastingsMover_fwd_hh
#define INCLUDED_protocols_canonical_sampling_SidechainMetropolisHastingsMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace protocols {
namespace canonical_sampling {

//Forwards and OP typedefs
class MetropolisHastingsMover;
typedef utility::pointer::shared_ptr< MetropolisHastingsMover > MetropolisHastingsMoverOP;
typedef utility::pointer::shared_ptr< MetropolisHastingsMover const > MetropolisHastingsMoverCOP;
typedef utility::pointer::weak_ptr< MetropolisHastingsMover > MetropolisHastingsMoverAP;
typedef utility::pointer::weak_ptr< MetropolisHastingsMover const > MetropolisHastingsMoverCAP;

} //moves
} //protocols

#endif //INCLUDED_protocols_canonical_sampling_MetropolisHastingsMover_FWD_HH
