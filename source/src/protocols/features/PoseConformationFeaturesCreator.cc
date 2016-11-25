// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/feature/PoseConformationFeaturesCreator.hh
/// @brief  Header for PoseConformationFeaturesCreator for the PoseConformationFeatures load-time factory registration scheme
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/PoseConformationFeaturesCreator.hh>

// Package Headers

#include <protocols/features/FeaturesReporterCreator.hh>

#include <protocols/features/PoseConformationFeatures.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace features {

// XRW TEMP PoseConformationFeaturesCreator::PoseConformationFeaturesCreator() {}
// XRW TEMP PoseConformationFeaturesCreator::~PoseConformationFeaturesCreator() = default;
// XRW TEMP FeaturesReporterOP PoseConformationFeaturesCreator::create_features_reporter() const {
// XRW TEMP  return FeaturesReporterOP( new PoseConformationFeatures );
// XRW TEMP }

// XRW TEMP std::string PoseConformationFeaturesCreator::type_name() const {
// XRW TEMP  return "PoseConformationFeatures";
// XRW TEMP }

} //namespace features
} //namespace protocols
