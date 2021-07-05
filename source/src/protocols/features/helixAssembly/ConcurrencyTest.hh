// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ConcurrencyTest.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_features_helixAssembly_ConcurrencyTest_hh
#define INCLUDED_protocols_features_helixAssembly_ConcurrencyTest_hh


#include <protocols/features/helixAssembly/ConcurrencyTest.fwd.hh>

//Core
#include <core/types.hh>

//Devel

//Utility and basic
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>

//C++
#include <string>

//External Headers

//Basic

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#include <protocols/features/FeaturesReporter.hh> // AUTO IWYU For FeaturesReporter

namespace protocols {
namespace features {
namespace helixAssembly {

class ConcurrencyTest : public protocols::features::FeaturesReporter
{

public:

	ConcurrencyTest(){}


	/// @brief generate the table schemas and write them to the database
	void
	write_schema_to_db(utility::sql_database::sessionOP db_session) const override;

	/// @brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1<bool> const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session
	) override;

	std::string
	type_name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

};

}
}
}
#endif
