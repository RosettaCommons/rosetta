// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/helical_bundle/BackboneGridSamplerHelper.fwd.hh
/// @brief  Owning pointers and whatnot for BackboneGridSamplerHelper class.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_protocols_helical_bundle_BackboneGridSamplerHelper_fwd_hh
#define INCLUDED_protocols_helical_bundle_BackboneGridSamplerHelper_fwd_hh

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>

// C++ headers

namespace protocols {
namespace helical_bundle {

class BackboneGridSamplerHelper;

typedef  utility::pointer::weak_ptr< BackboneGridSamplerHelper >  BackboneGridSamplerHelperAP;
typedef  utility::pointer::weak_ptr< BackboneGridSamplerHelper const >  BackboneGridSamplerHelperCAP;
typedef  utility::pointer::shared_ptr< BackboneGridSamplerHelper >  BackboneGridSamplerHelperOP;
typedef  utility::pointer::shared_ptr< BackboneGridSamplerHelper const >  BackboneGridSamplerHelperCOP;

typedef  utility::vector1< BackboneGridSamplerHelperOP >  BackboneGridSamplerHelperOPs;
typedef  utility::vector1< BackboneGridSamplerHelperCOP >  BackboneGridSamplerHelperCOPs;
typedef  utility::vector1< BackboneGridSamplerHelperCAP >  BackboneGridSamplerHelperCAPs;

} // namespace helical_bundle
} // namespace protocols

#endif // INCLUDED_protocols_helical_bundle_BackboneGridSamplerHelper_fwd_hh
