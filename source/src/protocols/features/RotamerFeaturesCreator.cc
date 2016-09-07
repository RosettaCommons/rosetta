// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/feature/RotamerFeaturesCreator.hh
/// @brief  Header for RotamerFeaturesCreator for the RotamerFeatures load-time factory registration scheme
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/RotamerFeaturesCreator.hh>

// Package Headers

#include <protocols/features/FeaturesReporterCreator.hh>

#include <protocols/features/RotamerFeatures.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace features {

RotamerFeaturesCreator::RotamerFeaturesCreator() {}
RotamerFeaturesCreator::~RotamerFeaturesCreator() = default;
FeaturesReporterOP RotamerFeaturesCreator::create_features_reporter() const {
	return FeaturesReporterOP( new RotamerFeatures );
}

std::string RotamerFeaturesCreator::type_name() const {
	return "RotamerFeatures";
}

} //namespace features
} //namespace protocols
