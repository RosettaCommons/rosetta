// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/helical_bundle/parameters/BundleParameters.fwd.hh
/// @brief  Owning pointers and whatnot for the class for holding parameters for parametric helical bundle generation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_protocols_helical_bundle_parameters_BundleParameters_fwd_hh
#define INCLUDED_protocols_helical_bundle_parameters_BundleParameters_fwd_hh

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>

// C++ headers

namespace protocols {
	namespace helical_bundle {
		namespace parameters {

			class BundleParameters;

			typedef  utility::pointer::weak_ptr< BundleParameters >  BundleParametersAP;
			typedef  utility::pointer::weak_ptr< BundleParameters const >  BundleParametersCAP;
			typedef  utility::pointer::shared_ptr< BundleParameters >  BundleParametersOP;
			typedef  utility::pointer::shared_ptr< BundleParameters const >  BundleParametersCOP;

			typedef  utility::vector1< BundleParametersOP >  BundleParametersOPs;
			typedef  utility::vector1< BundleParametersCOP >  BundleParametersCOPs;
			typedef  utility::vector1< BundleParametersCAP >  BundleParametersCAPs;

		} // namespace parameters
	} // namespace helical_bundle
} // namespace protocols

#endif // INCLUDED_protocols_helical_bundle_parameters_BundleParameters_fwd_hh
