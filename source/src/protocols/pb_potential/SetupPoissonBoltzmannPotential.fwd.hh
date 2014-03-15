// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/pb_potential/PoissonBoltzmannPotential.fwd.hh
/// @brief  setupPoissonBoltzmannPotential class forward delcaration
/// @author Sachko Honda (honda@apl.washington.edu)

#ifndef INCLUDED_protocols_pb_potential_SetupPoissonBoltzmannPotential_FWD_HH
#define INCLUDED_protocols_pb_potential_SetupPoissonBoltzmannPotential_FWD_HH

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace scoring {
class  SetupPoissonBoltzmannPotential;
typedef utility::pointer::owning_ptr< SetupPoissonBoltzmannPotential > SetupPoissonBoltzmannPotentialOP;
typedef utility::pointer::owning_ptr< SetupPoissonBoltzmannPotential const > SetupPoissonBoltzmannPotentialCOP;
typedef utility::pointer::access_ptr< SetupPoissonBoltzmannPotential > SetupPoissonBoltzmannPotentialAP;
typedef utility::pointer::access_ptr< SetupPoissonBoltzmannPotential const > SetupPoissonBoltzmannPotentialCAP;

}
}
#endif  // INCLUDED_protocols_pb_potential_SetupPoissonBoltzmannPotential_FWD_HH
