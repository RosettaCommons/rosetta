// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/energy_methods/ImplicitMembraneCoulomb.fwd.hh
/// @brief Minimal class for computing the depth- and membrane-dependent electrostatics energy
/// @author rfalford12 (rfalford12@gmail.com)

#ifndef INCLUDED_core_energy_methods_ImplicitMembraneCoulomb_fwd_hh
#define INCLUDED_core_energy_methods_ImplicitMembraneCoulomb_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace energy_methods {

class ImplicitMembraneCoulomb;

using ImplicitMembraneCoulombOP = utility::pointer::shared_ptr< ImplicitMembraneCoulomb >;
using ImplicitMembraneCoulombCOP = utility::pointer::shared_ptr< ImplicitMembraneCoulomb const >;

} //energy_methods
} //core

#endif //INCLUDED_core_energy_methods_ImplicitMembraneCoulomb_fwd_hh
