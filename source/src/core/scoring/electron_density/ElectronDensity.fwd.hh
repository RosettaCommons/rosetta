// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

//// @file   core/scoring/methods/electron_density/ElectronDensity.fwd.hh
/// @brief  Electron density forward declarations
/// @author Frank DiMaio

#ifndef INCLUDED_core_scoring_electron_density_ElectronDensity_fwd_hh
#define INCLUDED_core_scoring_electron_density_ElectronDensity_fwd_hh

//Utility headers
//too bad we can't use owning_ptr.fwd.hh
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace electron_density {

class ElectronDensity;
typedef utility::pointer::shared_ptr< ElectronDensity > ElectronDensityOP;
typedef utility::pointer::shared_ptr< ElectronDensity const > ElectronDensityCOP;

} // namespace electron_density
} //namespace scoring
} //namespace core

#endif
