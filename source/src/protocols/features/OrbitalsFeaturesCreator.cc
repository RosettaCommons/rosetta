// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/feature/OrbitalsFeaturesCreator.hh
/// @brief  Header for OrbitalsFeaturesCreator for the OrbitalsFeatures load-time factory registration scheme
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/OrbitalsFeaturesCreator.hh>

// Package Headers

#include <protocols/features/FeaturesReporterCreator.hh>

#include <protocols/features/OrbitalsFeatures.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace features {

OrbitalsFeaturesCreator::OrbitalsFeaturesCreator() {}
OrbitalsFeaturesCreator::~OrbitalsFeaturesCreator() {}
FeaturesReporterOP OrbitalsFeaturesCreator::create_features_reporter() const {
	return FeaturesReporterOP( new OrbitalsFeatures );
}

std::string OrbitalsFeaturesCreator::type_name() const {
	return "OrbitalsFeatures";
}

} //namespace features
} //namespace protocols
