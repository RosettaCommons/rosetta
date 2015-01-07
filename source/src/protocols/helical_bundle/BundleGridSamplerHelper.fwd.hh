// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/helical_bundle/BundleGridSamplerHelper.fwd.hh
/// @brief  Owning pointers and whatnot for BundleGridSamplerHelper class.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_protocols_helical_bundle_BundleGridSamplerHelper_fwd_hh
#define INCLUDED_protocols_helical_bundle_BundleGridSamplerHelper_fwd_hh

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>

// C++ headers

namespace protocols {
	namespace helical_bundle {

		class BundleGridSamplerHelper;

		typedef  utility::pointer::weak_ptr< BundleGridSamplerHelper >  BundleGridSamplerHelperAP;
		typedef  utility::pointer::weak_ptr< BundleGridSamplerHelper const >  BundleGridSamplerHelperCAP;
		typedef  utility::pointer::shared_ptr< BundleGridSamplerHelper >  BundleGridSamplerHelperOP;
		typedef  utility::pointer::shared_ptr< BundleGridSamplerHelper const >  BundleGridSamplerHelperCOP;

		typedef  utility::vector1< BundleGridSamplerHelperOP >  BundleGridSamplerHelperOPs;
		typedef  utility::vector1< BundleGridSamplerHelperCOP >  BundleGridSamplerHelperCOPs;
		typedef  utility::vector1< BundleGridSamplerHelperCAP >  BundleGridSamplerHelperCAPs;

	} // namespace helical_bundle
} // namespace protocols

#ifdef USEBOOSTSERIALIZE
#include <boost/serialization/access.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/map.hpp>
#endif

#endif // INCLUDED_protocols_helical_bundle_BundleGridSamplerHelper_fwd_hh
