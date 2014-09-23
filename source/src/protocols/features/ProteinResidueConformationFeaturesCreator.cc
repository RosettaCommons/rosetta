// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/feature/ProteinResidueConformationFeaturesCreator.hh
/// @brief  Header for ProteinResidueConformationFeaturesCreator for the ProteinResidueConformationFeatures load-time factory registration scheme
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/ProteinResidueConformationFeaturesCreator.hh>

// Package Headers

#include <protocols/features/FeaturesReporterCreator.hh>
// AUTO-REMOVED #include <protocols/features/FeaturesReporterFactory.hh>

#include <protocols/features/ProteinResidueConformationFeatures.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace features {

ProteinResidueConformationFeaturesCreator::ProteinResidueConformationFeaturesCreator() {}
ProteinResidueConformationFeaturesCreator::~ProteinResidueConformationFeaturesCreator() {}
FeaturesReporterOP ProteinResidueConformationFeaturesCreator::create_features_reporter() const {
	return FeaturesReporterOP( new ProteinResidueConformationFeatures );
}

std::string ProteinResidueConformationFeaturesCreator::type_name() const {
	return "ProteinResidueConformationFeatures";
}

} //namespace features
} //namespace protocols
