// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/features/FeaturesReporterFactory.hh
/// @brief
/// @author Matthew O'Meara (mattjomeara@gmail.com)


#ifndef INCLUDED_protocols_features_FeaturesReporterFactory_hh
#define INCLUDED_protocols_features_FeaturesReporterFactory_hh

// Unit Headers
#include <protocols/features/FeaturesReporterFactory.fwd.hh>

// Project Headers
#include <protocols/features/FeaturesReporter.hh>

// Platform Headers
#include <core/scoring/ScoreFunction.hh>

// C++ Headers
#include <string>

namespace protocols {
namespace features {

/// Create Features Reporters
class FeaturesReporterFactory {

	// Private constructor to make it singleton managed
	FeaturesReporterFactory();
	FeaturesReporterFactory(const FeaturesReporterFactory & src);

public:

	// Warning this is not called because of the singleton pattern
	virtual ~FeaturesReporterFactory();

	static
	FeaturesReporterOP
	create_features_reporter(
		std::string const & name,
		core::scoring::ScoreFunctionOP scfxn = NULL);

};

} // namespace
} // namespace

#endif
