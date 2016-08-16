// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/feature/ProteinBondGeometryFeaturesCreator.hh
/// @brief  Header for ProteinBondGeometryFeaturesCreator for the ProteinBondGeometryFeatures load-time factory registration scheme
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/ProteinBondGeometryFeaturesCreator.hh>

// Package Headers

#include <protocols/features/FeaturesReporterCreator.hh>

#include <protocols/features/ProteinBondGeometryFeatures.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace features {

ProteinBondGeometryFeaturesCreator::ProteinBondGeometryFeaturesCreator() {}
ProteinBondGeometryFeaturesCreator::~ProteinBondGeometryFeaturesCreator() {}
FeaturesReporterOP ProteinBondGeometryFeaturesCreator::create_features_reporter() const {
	return FeaturesReporterOP( new ProteinBondGeometryFeatures );
}

std::string ProteinBondGeometryFeaturesCreator::type_name() const {
	return "ProteinBondGeometryFeatures";
}

} //namespace features
} //namespace protocols
