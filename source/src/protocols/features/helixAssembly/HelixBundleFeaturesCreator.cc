// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file devel/helixAssembly/HelixBundleFeaturesCreator.cc
/// @brief
/// @author Tim Jacobs

// Unit Headers
#include <protocols/features/helixAssembly/HelixBundleFeaturesCreator.hh>

// Package Headers

#include <protocols/features/FeaturesReporterCreator.hh>

#include <protocols/features/helixAssembly/HelixBundleFeatures.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace features {
namespace helixAssembly {

HelixBundleFeaturesCreator::HelixBundleFeaturesCreator() {}
HelixBundleFeaturesCreator::~HelixBundleFeaturesCreator() {}
FeaturesReporterOP HelixBundleFeaturesCreator::create_features_reporter() const {
	return FeaturesReporterOP( new HelixBundleFeatures );
}

std::string HelixBundleFeaturesCreator::type_name() const {
	return "HelixBundleFeatures";
}

} //namespace helixAssembly
} //namespace features
} //namespace protocols
