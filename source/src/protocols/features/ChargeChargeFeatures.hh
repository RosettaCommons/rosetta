// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/ChargeChargeFeatures.hh
/// @brief  report comments stored with each pose
/// @author Joseph S Harrison

#ifndef INCLUDED_protocols_features_ChargeChargeFeatures_hh
#define INCLUDED_protocols_features_ChargeChargeFeatures_hh

// Unit Headers
#include <protocols/features/ChargeChargeFeatures.fwd.hh>
#include <protocols/features/FeaturesReporter.hh>

// Project Headers
#include <core/types.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.fwd.hh>


namespace protocols {
namespace features {

class ChargeChargeFeatures : public FeaturesReporter {
public:
	ChargeChargeFeatures();

	ChargeChargeFeatures(core::Length distance_cutoff);

	ChargeChargeFeatures(ChargeChargeFeatures const & src);

	~ChargeChargeFeatures() override= default;

	core::Length
	distance_cutoff() const { return distance_cutoff_; }

	void
	distance_cutoff(core::Length d) { distance_cutoff_ = d; }

	/// @brief return string with class name
	// XRW TEMP  std::string
	// XRW TEMP  type_name() const override;

	/// @brief generate the table schemas and write them to the database
	void
	write_schema_to_db(utility::sql_database::sessionOP db_session) const override;

	/// @brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const override;

	/// @brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const &,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session) override;

	std::string
	type_name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:

	core::Length distance_cutoff_;

};


} // features namespace
} // protocols namespace

#endif // include guard
