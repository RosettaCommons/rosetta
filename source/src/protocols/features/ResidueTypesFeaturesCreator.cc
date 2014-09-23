// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/feature/ResidueTypesFeaturesCreator.hh
/// @brief  Header for ResidueTypesFeaturesCreator for the ResidueTypesFeatures load-time factory registration scheme
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/ResidueTypesFeaturesCreator.hh>

// Package Headers

#include <protocols/features/FeaturesReporterCreator.hh>
// AUTO-REMOVED #include <protocols/features/FeaturesReporterFactory.hh>

#include <protocols/features/ResidueTypesFeatures.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace features {

ResidueTypesFeaturesCreator::ResidueTypesFeaturesCreator() {}
ResidueTypesFeaturesCreator::~ResidueTypesFeaturesCreator() {}
FeaturesReporterOP ResidueTypesFeaturesCreator::create_features_reporter() const {
	return FeaturesReporterOP( new ResidueTypesFeatures );
}

std::string ResidueTypesFeaturesCreator::type_name() const {
	return "ResidueTypesFeatures";
}

} //namespace features
} //namespace protocols
