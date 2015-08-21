// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/scoring/membrane/MPCbetaEnergy.fwd.hh
///
/// @brief  Membrane Environemnt CBeta Energy
/// @details One Body Term - Score packing density in the membrane. Scores centroids for within
///    6A and 12A radius. Derived from Membrane base potential and uses mpframework data
///    Last Modified: 4/2/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MPCbetaEnergy_fwd_hh
#define INCLUDED_core_scoring_membrane_MPCbetaEnergy_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace membrane {

/// @brief Membrane Potential Based CBeta Energy Term
class MPCbetaEnergy;
typedef utility::pointer::shared_ptr< MPCbetaEnergy > MPCbetaEnergyOP;
typedef utility::pointer::shared_ptr< MPCbetaEnergy const > MPCbetaEnergyCOP;

} // membrane
} // scoring
} // core


#endif // INCLUDED_core_scoring_membrane_MPCbetaEnergy_fwd_hh
