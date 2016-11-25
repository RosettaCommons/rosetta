// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/feature/GeometricSolvationFeaturesCreator.hh
/// @brief  Header for GeometricSolvationFeaturesCreator for the GeometricSolvationFeatures load-time factory registration scheme
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/GeometricSolvationFeaturesCreator.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// Package Headers

#include <protocols/features/FeaturesReporterCreator.hh>

#include <protocols/features/GeometricSolvationFeatures.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace features {

// XRW TEMP GeometricSolvationFeaturesCreator::GeometricSolvationFeaturesCreator() {}
// XRW TEMP GeometricSolvationFeaturesCreator::~GeometricSolvationFeaturesCreator() = default;
// XRW TEMP FeaturesReporterOP GeometricSolvationFeaturesCreator::create_features_reporter() const {
// XRW TEMP
// XRW TEMP  core::scoring::methods::EnergyMethodOptions options;
// XRW TEMP  return FeaturesReporterOP( new GeometricSolvationFeatures(options) );
// XRW TEMP }

// XRW TEMP std::string GeometricSolvationFeaturesCreator::type_name() const {
// XRW TEMP  return "GeometricSolvationFeatures";
// XRW TEMP }

} //namespace features
} //namespace protocols
