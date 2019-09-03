// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/helical_bundle_predict/HelicalBundlePredictApplication.fwd.hh
/// @brief The meat-and-potatoes for the helical_bundle_predict application, used to predict structures of helical bundles
/// made from canonical or noncanonical building-blocks.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_helical_bundle_predict_HelicalBundlePredictApplication_fwd_hh
#define INCLUDED_protocols_helical_bundle_predict_HelicalBundlePredictApplication_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace helical_bundle_predict {

class HelicalBundlePredictApplication;
class HelicalBundlePredictApplicationOptions;

typedef utility::pointer::shared_ptr< HelicalBundlePredictApplication > HelicalBundlePredictApplicationOP;
typedef utility::pointer::shared_ptr< HelicalBundlePredictApplication const > HelicalBundlePredictApplicationCOP;

typedef utility::pointer::shared_ptr< HelicalBundlePredictApplicationOptions > HelicalBundlePredictApplicationOptionsOP;
typedef utility::pointer::shared_ptr< HelicalBundlePredictApplicationOptions const > HelicalBundlePredictApplicationOptionsCOP;

} //protocols
} //helical_bundle_predict

#endif //INCLUDED_protocols_helical_bundle_predict_HelicalBundlePredictApplication_fwd_hh
