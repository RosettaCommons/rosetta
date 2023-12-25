// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/membrane/scoring/ElectricfieldLipidlayer.cc
/// @brief Implicit Lipid Membrane Model electrostatic energy due to the field created by lipid layers(one-body)
/// @author  Rituparna Samanta (rsamant2@jhu.edu)

#ifndef INCLUDED_protocols_membrane_scoring_ElectricfieldLipidlayer_fwd_hh
#define INCLUDED_protocols_membrane_scoring_ElectricfieldLipidlayer_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace membrane {
namespace scoring {

class ElectricfieldLipidlayer;
typedef utility::pointer::shared_ptr< ElectricfieldLipidlayer > ElectricfieldLipidlayerOP;
typedef utility::pointer::shared_ptr< ElectricfieldLipidlayer const > ElectricfieldLipidlayerCOP;

} // scoring
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_scoring_ElectricfieldLipidlayer_fwd_hh
