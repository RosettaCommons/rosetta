// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/FA_GrpElecEnergy.hh
/// @brief  Electrostatic energy with a distance-dependant dielectric
/// @author Hahnbeom Park


#ifndef INCLUDED_core_scoring_elec_FA_GrpElecEnergy_fwd_hh
#define INCLUDED_core_scoring_elec_FA_GrpElecEnergy_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace elec {

class FA_GrpElecEnergy;

typedef utility::pointer::shared_ptr< FA_GrpElecEnergy > FA_GrpElecEnergyOP;
typedef utility::pointer::shared_ptr< FA_GrpElecEnergy const > FA_GrpElecEnergyCOP;

class FAElecContextData;
typedef utility::pointer::shared_ptr< FAElecContextData > FAElecContextDataOP;
typedef utility::pointer::shared_ptr< FAElecContextData const > FAElecContextDataCOP;


}
}
}

#endif
