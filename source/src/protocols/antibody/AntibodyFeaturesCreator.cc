// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody_design/AntibodyFeaturesCreator.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit Headers
#include <protocols/antibody/AntibodyFeaturesCreator.hh>

// Package Headers

#include <protocols/features/FeaturesReporterCreator.hh>

#include <protocols/antibody/AntibodyFeatures.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace antibody {
using namespace protocols::features;

// XRW TEMP AntibodyFeaturesCreator::AntibodyFeaturesCreator() {}
// XRW TEMP AntibodyFeaturesCreator::~AntibodyFeaturesCreator() = default;
// XRW TEMP FeaturesReporterOP AntibodyFeaturesCreator::create_features_reporter() const {
// XRW TEMP  return FeaturesReporterOP( new AntibodyFeatures );
// XRW TEMP }

// XRW TEMP std::string AntibodyFeaturesCreator::type_name() const {
// XRW TEMP  return "AntibodyFeatures";
// XRW TEMP }

} //antibody
} //protocols
