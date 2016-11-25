// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/ScreeningFeatures.hh
/// @brief  report JSON object with information needed for vHTS postprocessing.
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_features_ScreeningFeatures_hh
#define INCLUDED_protocols_features_ScreeningFeatures_hh

#include <protocols/features/ScreeningFeatures.fwd.hh>
#include <protocols/features/FeaturesReporter.hh>

#include <utility/sql_database/DatabaseSessionManager.fwd.hh>
#include <utility/json_spirit/json_spirit_value.h>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.fwd.hh>


namespace protocols {
namespace features {

class ScreeningFeatures : public protocols::features::FeaturesReporter {

public:

	ScreeningFeatures();
	ScreeningFeatures(ScreeningFeatures const & src);

	~ScreeningFeatures() override;

	// XRW TEMP  std::string type_name() const override;

	void write_schema_to_db(utility::sql_database::sessionOP db_session) const override;

	utility::vector1<std::string> features_reporter_dependencies() const override;


	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & /*relevant_residues*/,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session) override;

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & /*data*/,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/) override;

	std::string
	type_name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	std::vector<utility::json_spirit::Pair> get_desriptor_data() const;

private:
	std::string chain_;
	utility::vector1<std::string> descriptors_;


};

}
}

#endif
