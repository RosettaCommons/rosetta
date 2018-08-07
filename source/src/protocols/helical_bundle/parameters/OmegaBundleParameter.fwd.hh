// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/helical_bundle/parameters/OmegaBundleParameter.fwd.hh
/// @brief Forward declarations for a class for the omega paremeter, derived from the generic RealValuedParameter class.
/// @details Omega has a few additional options that can be configured, for pitch-copying vs. value-copying.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_protocols_helical_bundle_parameters_OmegaBundleParameter_fwd_hh
#define INCLUDED_protocols_helical_bundle_parameters_OmegaBundleParameter_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Forward
namespace protocols {
namespace helical_bundle {
namespace parameters {

class OmegaBundleParameter;

typedef utility::pointer::shared_ptr< OmegaBundleParameter > OmegaBundleParameterOP;
typedef utility::pointer::shared_ptr< OmegaBundleParameter const > OmegaBundleParameterCOP;

} //protocols
} //helical_bundle
} //parameters

#endif //INCLUDED_protocols_helical_bundle_parameters_OmegaBundleParameter_fwd_hh
