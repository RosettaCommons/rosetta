// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/scoring/membrane/MPHelicalityEnergy.fwd.hh
///
/// @brief  Fullatom and centroid level smooth membrane non-helicality penalty
/// @details
/// @FlesihmanLab
/// @Last Modified: 20/2/17
///
/// @author Jonathan Weinstein (jonathan.weinstein@weizmann.ac.il)
/// @author Assaf Elazar
/// @author Sarel Fleishman

#ifndef INCLUDED_core_scoring_membrane_MPHelicalityEnergy_fwd_hh
#define INCLUDED_core_scoring_membrane_MPHelicalityEnergy_fwd_hh

// Utility Methods
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace membrane {

class MPHelicalityEnergy;
typedef utility::pointer::shared_ptr< MPHelicalityEnergy > MPHelicalityEnergyOP;
typedef utility::pointer::shared_ptr< MPHelicalityEnergy const > MPHelicalityEnergyCOP;

} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_membrane_MPHelicalityEnergy_fwd_hh
