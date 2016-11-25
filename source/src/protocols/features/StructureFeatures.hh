// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/StructureFeatures.hh
/// @brief  report structure identifiers to a features database
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_StructureFeatures_hh
#define INCLUDED_protocols_features_StructureFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/StructureFeatures.fwd.hh>

//External

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>
#include <utility/exit.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace features {

class StructureFeatures : public protocols::features::FeaturesReporter {
public:
	StructureFeatures();

	StructureFeatures( StructureFeatures const & src );

	~StructureFeatures() override;

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
		core::pose::Pose const &,
		utility::vector1< bool > const &,
		StructureID,
		utility::sql_database::sessionOP
	) override {
		utility_exit_with_message(
			"The Structure features reporter is a special feature reporter "
			" Please use the other interface for report_features");
		return 0;
	}


	/// @brief collect all the feature data for the pose use
	///This version allows the tag and the input tag to be specificed
	StructureID
	report_features(
		core::Size batch_id,
		utility::sql_database::sessionOP db_session,
		std::string const & tag,
		std::string const & input_tag
	);


	void mark_structure_as_sampled(
		core::Size batch_id,
		std::string const & tag,
		std::string const & input_tag,
		utility::sql_database::sessionOP db_session
	);

	void delete_record(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session
	) override;

	void
	load_into_pose(
		utility::sql_database::sessionOP db_session,
		StructureID struct_id,
		core::pose::Pose & pose) override;

	void
	load_tag(
		utility::sql_database::sessionOP db_session,
		StructureID struct_id,
		core::pose::Pose & pose);


	StructureID
	get_struct_id(
		utility::sql_database::sessionOP db_session,
		std::string const & tag,
		core::Size const & protocol_id
	);

	std::string
	type_name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


};

} // namespace
} // namespace

#endif // include guard
