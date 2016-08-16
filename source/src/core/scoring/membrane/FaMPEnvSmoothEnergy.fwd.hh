// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/scoring/membrane/FaMPEnvSmoothEnergy.fwd.hh
///
/// @brief  Fullatom Smoothed Membrane Environment Energy
/// @details Updated residue-environment energy (fullatom) by Vladmir in 2010 - smoothed
///    derivatives based on updated statistics. Adapted for mpframework by Rebecca
///    @GrayLab.
///    Last Modified: 7/6/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_FaMPEnvSmoothEnergy_fwd_hh
#define INCLUDED_core_scoring_membrane_FaMPEnvSmoothEnergy_fwd_hh

// Utility Methods
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace membrane {

class FaMPEnvSmoothEnergy;
typedef utility::pointer::shared_ptr< FaMPEnvSmoothEnergy > FaMPEnvSmoothEnergyOP;
typedef utility::pointer::shared_ptr< FaMPEnvSmoothEnergy const > FaMPEnvSmoothEnergyCOP;

} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_membrane_FaMPEnvSmoothEnergy_fwd_hh
