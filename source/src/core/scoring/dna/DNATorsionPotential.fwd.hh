// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/DNATorsionPotential.fwd.hh
/// @brief  DNATorsionPotential potential class forward delcaration
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_scoring_dna_DNATorsionPotential_FWD_HH
#define INCLUDED_core_scoring_dna_DNATorsionPotential_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace dna {

class DNATorsionPotential;
typedef utility::pointer::shared_ptr< DNATorsionPotential > DNATorsionPotentialOP;
typedef utility::pointer::shared_ptr< DNATorsionPotential const > DNATorsionPotentialCOP;

}
}
}

#endif
