// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file BatchFeaturesCreator.cc
///
/// @brief

/// @author Tim Jacobs

// Unit Headers
#include <protocols/features/BatchFeaturesCreator.hh>

// Package Headers

#include <protocols/features/FeaturesReporterCreator.hh>
// AUTO-REMOVED #include <protocols/features/FeaturesReporterFactory.hh>

#include <protocols/features/BatchFeatures.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace features {

		BatchFeaturesCreator::BatchFeaturesCreator() {}
		BatchFeaturesCreator::~BatchFeaturesCreator() {}
		FeaturesReporterOP BatchFeaturesCreator::create_features_reporter() const {
				return new BatchFeatures;
		}

		std::string BatchFeaturesCreator::type_name() const {
				return "BatchFeatures";
		}

} //namespace features
} //namespace protocols
