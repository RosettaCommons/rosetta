// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/feature/SaltBridgeFeaturesCreator.hh
/// @brief  Header for SaltBridgeFeaturesCreator for the SaltBridgeFeatures load-time factory registration scheme
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/SaltBridgeFeaturesCreator.hh>

// Package Headers

#include <protocols/features/FeaturesReporterCreator.hh>

#include <protocols/features/SaltBridgeFeatures.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace features {

// XRW TEMP SaltBridgeFeaturesCreator::SaltBridgeFeaturesCreator() {}
// XRW TEMP SaltBridgeFeaturesCreator::~SaltBridgeFeaturesCreator() = default;
// XRW TEMP FeaturesReporterOP SaltBridgeFeaturesCreator::create_features_reporter() const {
// XRW TEMP  return FeaturesReporterOP( new SaltBridgeFeatures );
// XRW TEMP }

// XRW TEMP std::string SaltBridgeFeaturesCreator::type_name() const {
// XRW TEMP  return "SaltBridgeFeatures";
// XRW TEMP }

} //namespace features
} //namespace protocols
