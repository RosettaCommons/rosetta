// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/feature/WaterFeaturesCreator.hh
/// @brief  Header for WaterFeaturesCreator for the WaterFeatures load-time factory registration scheme
/// @author Kevin Houlihan

// Unit Headers
#include <protocols/features/WaterFeaturesCreator.hh>

// Package Headers

#include <protocols/features/FeaturesReporterCreator.hh>
// AUTO-REMOVED #include <protocols/features/FeaturesReporterFactory.hh>

#include <protocols/features/WaterFeatures.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace features {

WaterFeaturesCreator::WaterFeaturesCreator() {}
WaterFeaturesCreator::~WaterFeaturesCreator() {}
FeaturesReporterOP WaterFeaturesCreator::create_features_reporter() const {
	return FeaturesReporterOP( new WaterFeatures );
}

std::string WaterFeaturesCreator::type_name() const {
	return "WaterFeatures";
}

} //namespace features
} //namespace protocols
