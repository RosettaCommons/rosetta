// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/beta_barrel/parameters/BarrelParameters.fwd.hh
/// @brief  Forward declarations for BarrelParameters.
/// @author Andy Watkins

#ifndef INCLUDED_protocols_beta_barrel_parameters_BarrelParameters_fwd_hh
#define INCLUDED_protocols_beta_barrel_parameters_BarrelParameters_fwd_hh

#include <utility/vector1.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace beta_barrel {
namespace parameters {

class BarrelParameters;

typedef utility::pointer::weak_ptr< BarrelParameters > BarrelParametersAP;
typedef utility::pointer::weak_ptr< BarrelParameters const > BarrelParametersCAP;
typedef utility::pointer::shared_ptr< BarrelParameters > BarrelParametersOP;
typedef utility::pointer::shared_ptr< BarrelParameters const > BarrelParametersCOP;

typedef utility::vector1< BarrelParametersOP > BarrelParametersOPs;
typedef utility::vector1< BarrelParametersCOP > BarrelParametersCOPs;

} // namespace parameters
} // namespace beta_barrel
} // namespace protocols

#endif
