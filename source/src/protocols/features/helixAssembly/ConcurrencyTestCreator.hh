// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ConcurrencyTestCreator.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_features_helixAssembly_ConcurrencyTestCreator_HH
#define INCLUDED_protocols_features_helixAssembly_ConcurrencyTestCreator_HH

// Unit Headers
#include <protocols/features/FeaturesReporterCreator.hh>

#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace features {
namespace helixAssembly {

/// @brief creator for the HelixBundleFeatures class
class ConcurrencyTestCreator : public FeaturesReporterCreator
{
public:
	// XRW TEMP  ConcurrencyTestCreator();
	// XRW TEMP  virtual ~ConcurrencyTestCreator();

	// XRW TEMP  virtual FeaturesReporterOP create_features_reporter() const;
	// XRW TEMP  virtual std::string type_name() const;
	protocols::features::FeaturesReporterOP create_features_reporter() const override;
	std::string type_name() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

} //namespace helixAssembly
} //namespace features
} //namespace protocols


#endif
