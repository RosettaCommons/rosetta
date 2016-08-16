// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/MMBondAngleEnergy.fwd.hh
/// @brief  molecular mechanics bond angle energy forward declaration
/// @author Colin A. Smith (colin.smith@ucsf.edu)

#ifndef INCLUDED_core_scoring_methods_MMBondAngleEnergy_fwd_hh
#define INCLUDED_core_scoring_methods_MMBondAngleEnergy_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace scoring {
namespace methods {

class MMBondAngleEnergy;

typedef  utility::pointer::weak_ptr< MMBondAngleEnergy > MMBondAngleEnergyAP;
typedef  utility::pointer::weak_ptr< MMBondAngleEnergy const > MMBondAngleEnergyCAP;
typedef  utility::pointer::shared_ptr< MMBondAngleEnergy > MMBondAngleEnergyOP;
typedef  utility::pointer::shared_ptr< MMBondAngleEnergy const > MMBondAngleEnergyCOP;

} // namespace methods
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_methods_MMBondAngleEnergy_HH
