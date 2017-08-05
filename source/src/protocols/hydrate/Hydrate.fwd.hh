// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file src/protocols/hydrate/Hydrate.cc
/// @brief The Hydrate Protocol
/// @detailed
/// @author Joaquin Ambia, Jason K. Lai

#ifndef INCLUDED_protocols_hydrate_Hydrate_fwd_hh
#define INCLUDED_protocols_hydrate_Hydrate_fwd_hh

// Utility header
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace hydrate {

/// @brief
class Hydrate;

//typedef utility::pointer::owning_ptr< Hydrate > HydrateOP;
//typedef utility::pointer::owning_ptr< Hydrate const > HydrateCOP;
typedef utility::pointer::shared_ptr< Hydrate > HydrateOP;
typedef utility::pointer::shared_ptr< Hydrate const > HydrateCOP;

} // namespace hydrate
} // namespace protocols

#endif  // INCLUDED_protocols_hydrate_Hydrate_FWD_HH

