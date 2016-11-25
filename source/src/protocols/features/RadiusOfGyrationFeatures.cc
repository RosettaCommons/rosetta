// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/RadiusOfGyrationFeatures.cc
/// @brief  report the radius of gyration to features Statistics Scientific Benchmark
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/features/RadiusOfGyrationFeatures.hh>
#include <core/scoring/methods/RG_Energy_Fast.hh>

// Project Headers
#include <core/types.hh>
#include <basic/database/sql_utils.hh>

// Utility Headers
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>

#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>

// External Headers
#include <cppdb/frontend.h>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/features/feature_schemas.hh>
#include <protocols/features/RadiusOfGyrationFeaturesCreator.hh>

namespace protocols {
namespace features {

using std::string;
using core::Size;
using core::pose::Pose;
using core::scoring::methods::RG_Energy_Fast;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;

RadiusOfGyrationFeatures::RadiusOfGyrationFeatures(){}

RadiusOfGyrationFeatures::RadiusOfGyrationFeatures( RadiusOfGyrationFeatures const & ) :
	FeaturesReporter()
{}

RadiusOfGyrationFeatures::~RadiusOfGyrationFeatures()= default;

// XRW TEMP string
// XRW TEMP RadiusOfGyrationFeatures::type_name() const { return "RadiusOfGyrationFeatures"; }

void
RadiusOfGyrationFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	write_radius_of_gyration_table_schema(db_session);
}

void
RadiusOfGyrationFeatures::write_radius_of_gyration_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column radius_of_gyration("radius_of_gyration", DbDataTypeOP( new DbReal() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(struct_id);
	vector1< std::string > reference_columns;
	reference_columns.push_back("struct_id");
	ForeignKey foreign_key(foreign_key_columns, "structures", reference_columns, true);

	Schema table("radius_of_gyration", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(radius_of_gyration);

	table.write(db_session);
}

utility::vector1<std::string>
RadiusOfGyrationFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("StructureFeatures");
	return dependencies;
}

Size
RadiusOfGyrationFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID struct_id,
	sessionOP db_session
){
	RG_Energy_Fast rg;

	std::string statement_string =  "INSERT INTO radius_of_gyration (struct_id, radius_of_gyration) VALUES (?,?);";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);
	stmt.bind(2,rg.calculate_rg_score(pose, relevant_residues));
	basic::database::safely_write_to_database(stmt);
	return 0;
}

std::string RadiusOfGyrationFeatures::type_name() const {
	return class_name();
}

std::string RadiusOfGyrationFeatures::class_name() {
	return "RadiusOfGyrationFeatures";
}

void RadiusOfGyrationFeatures::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::features::xsd_type_definition_w_attributes( xsd, class_name(), "Record radius of gyration for each structure as a feature", attlist );
}

std::string RadiusOfGyrationFeaturesCreator::type_name() const {
	return RadiusOfGyrationFeatures::class_name();
}

protocols::features::FeaturesReporterOP
RadiusOfGyrationFeaturesCreator::create_features_reporter() const {
	return protocols::features::FeaturesReporterOP( new RadiusOfGyrationFeatures );
}

void RadiusOfGyrationFeaturesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RadiusOfGyrationFeatures::provide_xml_schema( xsd );
}


} // features namesapce
} // protocols namespace
