// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/helical_bundle/parameters/BundleParametersSet.fwd.hh
/// @brief  Owning pointers and whatnot for the class for holding sets of parameters for helical bundle backbone generation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_protocols_helical_bundle_parameters_BundleParametersSet_fwd_hh
#define INCLUDED_protocols_helical_bundle_parameters_BundleParametersSet_fwd_hh

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>

// C++ headers

namespace protocols {
	namespace helical_bundle {
		namespace parameters {

			class BundleParametersSet;

			typedef  utility::pointer::weak_ptr< BundleParametersSet >  BundleParametersSetAP;
			typedef  utility::pointer::weak_ptr< BundleParametersSet const >  BundleParametersSetCAP;
			typedef  utility::pointer::shared_ptr< BundleParametersSet >  BundleParametersSetOP;
			typedef  utility::pointer::shared_ptr< BundleParametersSet const >  BundleParametersSetCOP;

			typedef  utility::vector1< BundleParametersSetOP >  BundleParametersSetOPs;
			typedef  utility::vector1< BundleParametersSetCOP >  BundleParametersSetCOPs;
			typedef  utility::vector1< BundleParametersSetCAP >  BundleParametersSetCAPs;

		} // namespace parameters
	} // namespace helical_bundle
} // namespace protocols

#endif // INCLUDED_protocols_helical_bundle_parameters_BundleParametersSet_fwd_hh
