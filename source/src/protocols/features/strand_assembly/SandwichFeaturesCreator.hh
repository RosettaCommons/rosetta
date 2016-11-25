// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/features/strand_assembly/SandwichFeaturesCreator.hh
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_features_strand_assembly_SandwichFeaturesCreator_hh
#define INCLUDED_protocols_features_strand_assembly_SandwichFeaturesCreator_hh

// Unit Headers
#include <protocols/features/FeaturesReporterCreator.hh>

#include <core/types.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace features {
namespace strand_assembly {

/// @brief creator for the SandwichFeatures class
class SandwichFeaturesCreator : public FeaturesReporterCreator
{
public:
	// XRW TEMP  SandwichFeaturesCreator();
	// XRW TEMP  virtual ~SandwichFeaturesCreator();

	// XRW TEMP  virtual FeaturesReporterOP create_features_reporter() const;
	// XRW TEMP  virtual std::string type_name() const;
	protocols::features::FeaturesReporterOP create_features_reporter() const override;
	std::string type_name() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

} //namespace strand_assembly
} //namespace features
} //namespace protocols

#endif
