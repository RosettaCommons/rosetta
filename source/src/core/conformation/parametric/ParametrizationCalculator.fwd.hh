// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/parametric/ParametrizationCalculator.fwd.hh
/// @brief  Forward declarations for the ParametrizationCalculator base class, from which classes
/// that calculate particular parametrizations will be derived.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_core_conformation_parametric_ParametrizationCalculator_fwd_hh
#define INCLUDED_core_conformation_parametric_ParametrizationCalculator_fwd_hh

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>

// C++ headers

namespace core {
namespace conformation {
namespace parametric {

class ParametrizationCalculator;

typedef  utility::pointer::weak_ptr< ParametrizationCalculator >  ParametrizationCalculatorAP;
typedef  utility::pointer::weak_ptr< ParametrizationCalculator const >  ParametrizationCalculatorCAP;
typedef  utility::pointer::shared_ptr< ParametrizationCalculator >  ParametrizationCalculatorOP;
typedef  utility::pointer::shared_ptr< ParametrizationCalculator const >  ParametrizationCalculatorCOP;

typedef  utility::vector1< ParametrizationCalculatorOP >  ParametrizationCalculatorOPs;
typedef  utility::vector1< ParametrizationCalculatorCOP >  ParametrizationCalculatorCOPs;
typedef  utility::vector1< ParametrizationCalculatorCAP >  ParametrizationCalculatorCAPs;

} // namespace parametric
} // namespace conformation
} // namespace core

#endif // INCLUDED_core_conformation_parametric_ParametrizationCalculator_fwd_hh
