// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		core/scoring/membrane/MPLipoEnergy.fwd.hh
///
///	@brief		Membrane Lipophibicity Term
///	@details	Whole Structure Energy - Evaluate structure based on derived
///				lipophobicities from input in lips file.
///				Last Modified: 3/28/14
///
///	@author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MPLipoEnergy_fwd_hh
#define INCLUDED_core_scoring_membrane_MPLipoEnergy_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh> 

namespace core {
namespace scoring {
namespace membrane {
	
class MPLipoEnergy;
typedef utility::pointer::owning_ptr< MPLipoEnergy > MPLipoEnergyOP;
typedef utility::pointer::owning_ptr< MPLipoEnergy const > MPLipoEnergyCOP;
	
} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_membrane_MPLipoEnergy_fwd_hh
