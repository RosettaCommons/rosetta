// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/beta_barrel/BarrelParametrizationCalculator.fwd.hh
/// @brief  Forward declarations for BarrelParametrizationCalculator.
/// @author Andy Watkins

#ifndef INCLUDED_protocols_beta_barrel_BarrelParametrizationCalculator_fwd_hh
#define INCLUDED_protocols_beta_barrel_BarrelParametrizationCalculator_fwd_hh

#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace beta_barrel {

class BarrelParametrizationCalculator;

typedef utility::pointer::weak_ptr< BarrelParametrizationCalculator > BarrelParametrizationCalculatorAP;
typedef utility::pointer::weak_ptr< BarrelParametrizationCalculator const > BarrelParametrizationCalculatorCAP;
typedef utility::pointer::shared_ptr< BarrelParametrizationCalculator > BarrelParametrizationCalculatorOP;
typedef utility::pointer::shared_ptr< BarrelParametrizationCalculator const > BarrelParametrizationCalculatorCOP;

} // namespace beta_barrel
} // namespace protocols

#endif
