// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/feature/AtomTypesFeaturesCreator.hh
/// @brief  Header for AtomTypesFeaturesCreator for the AtomTypesFeatures load-time factory registration scheme
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/AtomTypesFeaturesCreator.hh>

// Package Headers

#include <protocols/features/FeaturesReporterCreator.hh>
// AUTO-REMOVED #include <protocols/features/FeaturesReporterFactory.hh>

#include <protocols/features/AtomTypesFeatures.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace features {

AtomTypesFeaturesCreator::AtomTypesFeaturesCreator() {}
AtomTypesFeaturesCreator::~AtomTypesFeaturesCreator() {}
FeaturesReporterOP AtomTypesFeaturesCreator::create_features_reporter() const {
	return FeaturesReporterOP( new AtomTypesFeatures );
}

std::string AtomTypesFeaturesCreator::type_name() const {
	return "AtomTypesFeatures";
}

} //namespace features
} //namespace protocols
