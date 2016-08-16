// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/feature/BetaTurnDetectionFeaturesCreator.hh
/// @brief  Header for BetaTurnDetectionFeaturesCreator for the BetaTurnDetectionFeatures load-time factory registration scheme
/// @author Brian Weitzner

// Unit Headers
#include <protocols/features/BetaTurnDetectionFeaturesCreator.hh>

// Package Headers

#include <protocols/features/FeaturesReporterCreator.hh>

#include <protocols/features/BetaTurnDetectionFeatures.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace features {

BetaTurnDetectionFeaturesCreator::BetaTurnDetectionFeaturesCreator() {}
BetaTurnDetectionFeaturesCreator::~BetaTurnDetectionFeaturesCreator() {}
FeaturesReporterOP BetaTurnDetectionFeaturesCreator::create_features_reporter() const {
	return FeaturesReporterOP( new BetaTurnDetectionFeatures );
}

std::string BetaTurnDetectionFeaturesCreator::type_name() const {
	return "BetaTurnDetectionFeatures";
}

} //namespace features
} //namespace protocols
