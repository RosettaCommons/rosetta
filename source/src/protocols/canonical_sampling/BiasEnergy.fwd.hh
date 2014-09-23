// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
//  vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   /protocols/canonical_sampling/BiasEnergy.fwd.hh
/// @brief
/// @author Oliver Lange  ( oliver.lange@tum.de )

#ifndef INCLUDED_protocols_canonical_sampling_BiasEnergy_fwd_hh
#define INCLUDED_protocols_canonical_sampling_BiasEnergy_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace protocols {
namespace canonical_sampling {

// Forward
class BiasEnergy;

typedef utility::pointer::shared_ptr< BiasEnergy > BiasEnergyOP;
typedef utility::pointer::shared_ptr< BiasEnergy const > BiasEnergyCOP;
typedef utility::pointer::weak_ptr< BiasEnergy > BiasEnergyAP;
typedef utility::pointer::weak_ptr< BiasEnergy const > BiasEnergyCAP;

// Forward


} // namespace canonical_sampling
} // namespace protocols


#endif // INCLUDED_protocols_canonical_sampling_BiasEnergy_FWD_HH

