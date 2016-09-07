// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/feature/ScreeningFeaturesCreator.cc
/// @brief
/// @author Sam DeLuca

// Unit Headers
#include <protocols/features/ScreeningFeaturesCreator.hh>

// Package Headers
#include <protocols/features/FeaturesReporterCreator.hh>

#include <protocols/features/ScreeningFeatures.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace features {

ScreeningFeaturesCreator::ScreeningFeaturesCreator() {}
ScreeningFeaturesCreator::~ScreeningFeaturesCreator() = default;
FeaturesReporterOP ScreeningFeaturesCreator::create_features_reporter() const {
	return FeaturesReporterOP( new ScreeningFeatures );
}

std::string ScreeningFeaturesCreator::type_name() const {
	return "ScreeningFeatures";
}

} //namespace features
} //namespace protocols
