// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/ScoreTypeFeatures.cc
/// @brief  report protocol level features to features statistics scientific benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/ScoreTypeFeatures.hh>
#include <protocols/features/util.hh>

//External

//Basic Headers
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>

// Platform Headers
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>
// #include <basic/database/sql_utils.hh> 09_11_2013, commented by Doonam due to double declaration

// External Headers

// C++ Headers
#include <sstream>
#include <string>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/features/feature_schemas.hh>
#include <protocols/features/ScoreTypeFeaturesCreator.hh>

#include <utility/string_util.hh> // AUTO IWYU For to_string

namespace protocols {
namespace features {

using std::string;
using std::stringstream;
using core::pose::Pose;
using core::scoring::EnergyMap;
using core::scoring::get_score_function;
using core::scoring::n_score_types;
using core::scoring::ScoreTypeManager;
using core::scoring::ScoreType;
using core::scoring::ScoreFunctionOP;
using core::Size;
using core::Real;
using utility::vector1;
using utility::sql_database::sessionOP;

ScoreTypeFeatures::ScoreTypeFeatures() :
	scfxn_(get_score_function())
{}

ScoreTypeFeatures::ScoreTypeFeatures(
	ScoreFunctionOP scfxn
) :
	scfxn_(std::move(scfxn))
{
	if ( scfxn_ == nullptr ) {
		utility_exit_with_message( "ScoreTypeFeatures may not be constructed with a null-pointer ScoreFunctionOP" );
	}
}

void
ScoreTypeFeatures::write_schema_to_db(
	sessionOP db_session
) const {

	using namespace basic::database::schema_generator;

	//******score_types******//
	Column batch_id("batch_id", utility::pointer::make_shared< DbInteger >(), false);
	Column score_type_id("score_type_id", utility::pointer::make_shared< DbInteger >(), false);
	Column score_type_name("score_type_name", utility::pointer::make_shared< DbText >(), false);

	utility::vector1<Column> pkey_cols;
	pkey_cols.push_back(batch_id);
	pkey_cols.push_back(score_type_id);

	Schema score_types("score_types", PrimaryKey(pkey_cols));
	score_types.add_column(batch_id);
	score_types.add_column(score_type_id);
	score_types.add_column(score_type_name);

	score_types.add_foreign_key(ForeignKey(batch_id, "batches", "batch_id", true));

	score_types.write(db_session);
}

utility::vector1<std::string>
ScoreTypeFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ProtocolFeatures");
	return dependencies;
}

Size
ScoreTypeFeatures::report_features(
	Pose const &,
	vector1< bool > const &,
	StructureID const struct_id,
	sessionOP db_session
){
	core::Size const batch_id(get_batch_id(struct_id, db_session));

	insert_score_type_rows(batch_id, db_session);
	return 0;
}


Size
ScoreTypeFeatures::report_features(
	core::Size const batch_id,
	sessionOP db_session
){
	insert_score_type_rows(batch_id, db_session);
	return 0;
}


void ScoreTypeFeatures::delete_record(
	StructureID,
	utility::sql_database::sessionOP ){}

void
ScoreTypeFeatures::insert_score_type_rows(
	core::Size protocol_id,
	sessionOP db_session
) {
	switch(db_session->get_db_mode()){
	case utility::sql_database::DatabaseMode::sqlite3:
	case utility::sql_database::DatabaseMode::mysql:
	case utility::sql_database::DatabaseMode::postgres : {
		using namespace basic::database;

		string const table_name("score_types");

		std::vector<std::string> const column_names{ "batch_id", "score_type_id", "score_type_name" };

		std::string protocol_id_s = utility::to_string (protocol_id);

		for ( core::Size score_type_id=1; score_type_id <= n_score_types; ++score_type_id ) {

			std::string score_type_id_s = utility::to_string (score_type_id);

			auto type(static_cast<ScoreType>(score_type_id));
			string const score_type( ScoreTypeManager::name_from_score_type(type) );

			std::vector<std::string> values;
			values.push_back(protocol_id_s);
			values.push_back(score_type_id_s);
			values.push_back(score_type);

			insert_or_ignore(table_name, column_names, values, db_session);
		}
	}
		break;

	default :
		utility_exit_with_message(
			"Unrecognized database mode: '" +
			name_from_database_mode(db_session->get_db_mode()) + "'");
		break;
	}
}

std::string ScoreTypeFeatures::type_name() const {
	return class_name();
}

std::string ScoreTypeFeatures::class_name() {
	return "ScoreTypeFeatures";
}

void ScoreTypeFeatures::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::features::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Report scoring type features",
		attlist);
}

std::string ScoreTypeFeaturesCreator::type_name() const {
	return ScoreTypeFeatures::class_name();
}

protocols::features::FeaturesReporterOP
ScoreTypeFeaturesCreator::create_features_reporter() const {
	return utility::pointer::make_shared< ScoreTypeFeatures >();
}

void ScoreTypeFeaturesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ScoreTypeFeatures::provide_xml_schema( xsd );
}


} // namesapce
} // namespace
