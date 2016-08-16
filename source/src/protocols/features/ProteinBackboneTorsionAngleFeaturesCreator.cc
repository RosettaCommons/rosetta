// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/feature/ProteinBackboneTorsionAngleFeaturesCreator.hh
/// @brief  Header for ProteinBackboneTorsionAngleFeaturesCreator for the ProteinBackboneTorsionAngleFeatures load-time factory registration scheme
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/ProteinBackboneTorsionAngleFeaturesCreator.hh>

// Package Headers

#include <protocols/features/FeaturesReporterCreator.hh>

#include <protocols/features/ProteinBackboneTorsionAngleFeatures.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace features {

ProteinBackboneTorsionAngleFeaturesCreator::ProteinBackboneTorsionAngleFeaturesCreator() {}
ProteinBackboneTorsionAngleFeaturesCreator::~ProteinBackboneTorsionAngleFeaturesCreator() {}
FeaturesReporterOP ProteinBackboneTorsionAngleFeaturesCreator::create_features_reporter() const {
	return FeaturesReporterOP( new ProteinBackboneTorsionAngleFeatures );
}

std::string ProteinBackboneTorsionAngleFeaturesCreator::type_name() const {
	return "ProteinBackboneTorsionAngleFeatures";
}

} //namespace features
} //namespace protocols
