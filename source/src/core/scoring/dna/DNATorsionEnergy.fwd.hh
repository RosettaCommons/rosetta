// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/dna/DNATorsionEnergy.fwd.hh
/// @brief
/// @author Jim Havranek


#ifndef INCLUDED_core_scoring_dna_DNATorsionEnergy_FWD_HH
#define INCLUDED_core_scoring_dna_DNATorsionEnergy_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace dna {

class DNATorsionEnergy;

typedef utility::pointer::shared_ptr< DNATorsionEnergy > DNATorsionEnergyOP;
typedef utility::pointer::shared_ptr< DNATorsionEnergy const > DNATorsionEnergyCOP;

} // dna
} // scoring
} // core


#endif
