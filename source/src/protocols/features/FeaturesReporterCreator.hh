// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/feature/FeaturesReporterCreator.hh
/// @brief  Base class for FeaturesReporterCreators for the FeaturesReporter load-time factory registration scheme
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_FeaturesReporterCreator_hh
#define INCLUDED_protocols_features_FeaturesReporterCreator_hh

// Unit Headers
#include <protocols/features/FeaturesReporterCreator.fwd.hh>

// Package Headers
#include <protocols/features/FeaturesReporter.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// c++ headers
#include <string>

namespace protocols {
namespace features {

/// @brief The Creator class is responsible for creating a particular
/// mover class.
class FeaturesReporterCreator : public utility::pointer::ReferenceCount
{
public:
	FeaturesReporterCreator() {}
	~FeaturesReporterCreator() override = default;

	virtual FeaturesReporterOP create_features_reporter() const = 0;
	virtual std::string type_name() const = 0;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const = 0;
};

} //namespace
} //namespace

#endif
