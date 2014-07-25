// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		core/scoring/membrane/MPPairEnergy.fwd.hh
///
///	@brief		Membrane Residue Pair Energy Term
///	@details	Two Body Term - score residue-residue interactions in the membrane. Derived from Membrane
///				base potential and uses mpframework data
///				Last Modified: 3/28/14
///
///	@author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MPPairEnergy_fwd_hh
#define INCLUDED_core_scoring_membrane_MPPairEnergy_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace membrane {
	
/// @brief Membrane Residue Pair Energy Term
class MPPairEnergy;
typedef utility::pointer::owning_ptr< MPPairEnergy > MPPairEnergyOP;
typedef utility::pointer::owning_ptr< MPPairEnergy const > MPPairEnergyCOP;

} // membrane
} // scoring
} // core


#endif // INCLUDED_core_scoring_membrane_MPPairEnergy_fwd_hh
