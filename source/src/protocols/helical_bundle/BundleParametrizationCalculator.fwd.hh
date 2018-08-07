// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/helical_bundle/BundleParametrizationCalculator.fwd.hh
/// @brief  Forward declarations for the BundleParametrizationCalculator class.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_protocols_helical_bundle_BundleParametrizationCalculator_fwd_hh
#define INCLUDED_protocols_helical_bundle_BundleParametrizationCalculator_fwd_hh

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>

// C++ headers

namespace protocols {
namespace helical_bundle {

class BundleParametrizationCalculator;

typedef  utility::pointer::weak_ptr< BundleParametrizationCalculator >  BundleParametrizationCalculatorAP;
typedef  utility::pointer::weak_ptr< BundleParametrizationCalculator const >  BundleParametrizationCalculatorCAP;
typedef  utility::pointer::shared_ptr< BundleParametrizationCalculator >  BundleParametrizationCalculatorOP;
typedef  utility::pointer::shared_ptr< BundleParametrizationCalculator const >  BundleParametrizationCalculatorCOP;

typedef  utility::vector1< BundleParametrizationCalculatorOP >  BundleParametrizationCalculatorOPs;
typedef  utility::vector1< BundleParametrizationCalculatorCOP >  BundleParametrizationCalculatorCOPs;
typedef  utility::vector1< BundleParametrizationCalculatorCAP >  BundleParametrizationCalculatorCAPs;

} // namespace helical_bundle
} // namespace protocols

#endif // INCLUDED_protocols_helical_bundle_BundleParametrizationCalculator_fwd_hh
