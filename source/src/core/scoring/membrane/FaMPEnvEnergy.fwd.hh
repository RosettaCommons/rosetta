// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/scoring/membrane/FaMPEnvEnergy.fwd.hh
///
/// @brief  LK-Type Membrane Environment Energy
/// @details Last Modified: 5/13/14
///
/// @author  Patrick Barth (Original)
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_FaMPEnvEnergy_fwd_hh
#define INCLUDED_core_scoring_membrane_FaMPEnvEnergy_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace membrane {

class FaMPEnvEnergy;
typedef utility::pointer::shared_ptr< FaMPEnvEnergy > FaMPEnvEnergyOP;
typedef utility::pointer::shared_ptr< FaMPEnvEnergy const > FaMPEnvEnergyCOP;

} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_methods_FaMPEnvEnergy_fwd_hh
