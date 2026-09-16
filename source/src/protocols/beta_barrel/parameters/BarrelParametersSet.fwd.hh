// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/beta_barrel/parameters/BarrelParametersSet.fwd.hh
/// @brief  Forward declarations for BarrelParametersSet.
/// @author Andy Watkins

#ifndef INCLUDED_protocols_beta_barrel_parameters_BarrelParametersSet_fwd_hh
#define INCLUDED_protocols_beta_barrel_parameters_BarrelParametersSet_fwd_hh

#include <utility/vector1.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace beta_barrel {
namespace parameters {

class BarrelParametersSet;

typedef utility::pointer::weak_ptr< BarrelParametersSet > BarrelParametersSetAP;
typedef utility::pointer::weak_ptr< BarrelParametersSet const > BarrelParametersSetCAP;
typedef utility::pointer::shared_ptr< BarrelParametersSet > BarrelParametersSetOP;
typedef utility::pointer::shared_ptr< BarrelParametersSet const > BarrelParametersSetCOP;

} // namespace parameters
} // namespace beta_barrel
} // namespace protocols

#endif
