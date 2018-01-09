// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/scoring/membrane/MPResidueLipophilicityEnergy.fwd.hh
///
/// @brief  Fullatom Smoothed Membrane Environment Energy
/// @details residue speicific enrgy by membrane depth, according to the Elazar
/// hydrophobicity scale
///    @FlesihmanLab.
///    Last Modified: 4/4/16
///
/// @author Jonathan Weinstein (jonathan.weinstein@weizmann.ac.il)
/// @author Assaf Elazar
/// @author Sarel Fleishman

#ifndef INCLUDED_core_scoring_membrane_MPResidueLipophilicityEnergy_fwd_hh
#define INCLUDED_core_scoring_membrane_MPResidueLipophilicityEnergy_fwd_hh

// Utility Methods
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace membrane {

class MPResidueLipophilicityEnergy;
typedef utility::pointer::shared_ptr< MPResidueLipophilicityEnergy > MPResidueLipophilicityEnergyOP;
typedef utility::pointer::shared_ptr< MPResidueLipophilicityEnergy const > MPResidueLipophilicityEnergyCOP;

} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_membrane_MPResidueLipophilicityEnergy_fwd_hh
