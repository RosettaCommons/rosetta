// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/MMBondLengthEnergy.fwd.hh
/// @brief  Molecular mechanics bond length score class
/// @author Frank DiMaio (based on Colin Smith's MMBondAngle potential)

#ifndef INCLUDED_core_scoring_methods_MMBondLengthEnergy_fwd_hh
#define INCLUDED_core_scoring_methods_MMBondLengthEnergy_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace scoring {
namespace methods {

class MMBondLengthEnergy;

typedef  utility::pointer::weak_ptr< MMBondLengthEnergy > MMBondLengthEnergyAP;
typedef  utility::pointer::weak_ptr< MMBondLengthEnergy const > MMBondLengthEnergyCAP;
typedef  utility::pointer::shared_ptr< MMBondLengthEnergy > MMBondLengthEnergyOP;
typedef  utility::pointer::shared_ptr< MMBondLengthEnergy const > MMBondLengthEnergyCOP;

} // namespace methods
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_methods_MMBondLengthEnergy_HH
