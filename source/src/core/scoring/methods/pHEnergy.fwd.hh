// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file   core/scoring/methods/pHEnergy.fwd.hh
/// @brief  pH energy forward declaration
/// @author krishna


#ifndef INCLUDED_core_scoring_methods_pHEnergy_fwd_hh
#define INCLUDED_core_scoring_methods_pHEnergy_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace methods {

class pHEnergy;

typedef utility::pointer::shared_ptr< pHEnergy > pHEnergyOP;

} // methods
} // scoring
} // core


#endif // INCLUDED_core_scoring_methods_pHEnergy_FWD_HH
