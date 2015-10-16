// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ScoreFunctionFeatures.cc
/// @brief  report the score parameters to a features database
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/ScoreFunctionFeatures.hh>
#include <protocols/features/util.hh>

//External

//Basic Headers
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/Constraint.hh>

// Platform Headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>
#include <basic/database/sql_utils.hh>
#include <utility/tag/Tag.hh>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <sstream>
#include <utility/excn/Exceptions.hh>
#include <string>


namespace protocols {
namespace features {

using std::string;
using std::stringstream;
using std::endl;
using basic::database::safely_prepare_statement;
using core::pose::Pose;
using core::graph::Graph;
using core::scoring::Energies;
using core::scoring::EnergyGraph;
using core::scoring::EnergyMap;
using core::scoring::EnergyEdge;
using core::scoring::get_score_function;
using core::scoring::get_score_functionName;
using core::scoring::hbond_sr_bb;
using core::scoring::hbond_lr_bb;
using core::scoring::n_score_types;
using core::scoring::ScoreFunction;
using core::scoring::ScoreFunctionOP;
using core::scoring::ScoreType;
using core::scoring::ScoreTypeManager;
using core::scoring::methods::EnergyMethodOptions;
using core::scoring::hbonds::HBondSet;
using core::scoring::hbonds::get_hbond_energies;
using core::scoring::EnergiesCacheableDataType::HBOND_SET;
using core::Size;
using core::Real;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;

ScoreFunctionFeatures::ScoreFunctionFeatures() :
	scfxn_(get_score_function()),
	scfxn_name_(get_score_functionName())
{}

ScoreFunctionFeatures::ScoreFunctionFeatures(
	ScoreFunctionOP scfxn,
	std::string const & scfxn_name
) :
	scfxn_(scfxn),
	scfxn_name_(scfxn_name)
{
	if ( scfxn_ == 0 ) {
		utility_exit_with_message( "ScoreFunctionFeatures may not be constructed with a null-pointer ScoreFunctionOP" );
	}
}

ScoreFunctionFeatures::ScoreFunctionFeatures(
	ScoreFunctionFeatures const & src
) :
	FeaturesReporter(),
	scfxn_(src.scfxn_),
	scfxn_name_(src.scfxn_name_)
{}

ScoreFunctionFeatures::~ScoreFunctionFeatures() {}

string
ScoreFunctionFeatures::type_name() const { return "ScoreFunctionFeatures"; }

void
ScoreFunctionFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	write_score_function_weights_table_schema(db_session);
	EnergyMethodOptions::write_score_function_method_options_table_schema(
		db_session);
}

void
ScoreFunctionFeatures::write_score_function_weights_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column batch_id("batch_id", DbDataTypeOP( new DbInteger() ), true);
	Column score_function_name("score_function_name", DbDataTypeOP( new DbText(255) ), true);
	Column score_type_id("score_type_id", DbDataTypeOP( new DbInteger() ), true);
	Column weight("weight", DbDataTypeOP( new DbReal() ), true);

	utility::vector1<Column> pkey_cols;
	pkey_cols.push_back(batch_id);
	pkey_cols.push_back(score_function_name);
	pkey_cols.push_back(score_type_id);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(batch_id);
	foreign_key_columns.push_back(score_type_id);
	vector1< string > reference_columns;
	reference_columns.push_back("batch_id");
	reference_columns.push_back("score_type_id");
	ForeignKey foreign_key(foreign_key_columns, "score_types", reference_columns, true);

	Schema table("score_function_weights", PrimaryKey(pkey_cols));
	table.add_foreign_key(foreign_key);
	table.add_column(weight);

	table.write(db_session);

}

vector1<string>
ScoreFunctionFeatures::features_reporter_dependencies() const {
	vector1<string> dependencies;
	dependencies.push_back("ScoreTypeFeatures");
	return dependencies;
}


void
ScoreFunctionFeatures::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	Pose const & /*pose*/
) {
	if ( tag->hasOption("scorefxn") ) {
		scfxn_name_ = tag->getOption<string>("scorefxn");
		scfxn_ = data.get_ptr<ScoreFunction>("scorefxns", scfxn_name_);
	}
}

Size
ScoreFunctionFeatures::report_features(
	Pose const &,
	vector1< bool > const &,
	StructureID const struct_id,
	sessionOP db_session
){
	Size const batch_id(get_batch_id(struct_id, db_session));

	insert_score_function_weights_rows(batch_id, db_session);
	scfxn_->energy_method_options().insert_score_function_method_options_rows(
		batch_id, scfxn_name_, db_session);
	return 0;

}

void
ScoreFunctionFeatures::delete_record(
	StructureID,
	utility::sql_database::sessionOP ){}

void
ScoreFunctionFeatures::insert_score_function_weights_rows(
	core::Size batch_id,
	sessionOP db_session
) const {

	string statement_string;
	switch(db_session->get_db_mode()){
	case utility::sql_database::DatabaseMode::sqlite3 :
		statement_string = "INSERT OR IGNORE INTO score_function_weights (batch_id, score_function_name, score_type_id, weight) VALUES (?,?,?,?);";
		break;
	case utility::sql_database::DatabaseMode::mysql:
	case utility::sql_database::DatabaseMode::postgres :
		statement_string = "INSERT IGNORE INTO score_function_weights (batch_id, score_function_name, score_type_id, weight) VALUES (?,?,?,?);";
		break;
	default :
		utility_exit_with_message(
			"Unrecognized database mode: '" +
			name_from_database_mode(db_session->get_db_mode()) + "'");
	}

	statement stmt(safely_prepare_statement(statement_string, db_session));

	for ( Size score_type_id=1; score_type_id <= n_score_types; ++score_type_id ) {
		ScoreType type(static_cast<ScoreType>(score_type_id));

		Real const weight( scfxn_->weights()[type] );
		if ( !weight ) continue;

		string const score_type( ScoreTypeManager::name_from_score_type(type) );
		stmt.bind(1, batch_id);
		stmt.bind(2, scfxn_name_);
		stmt.bind(3, score_type_id);
		stmt.bind(4, weight);
		basic::database::safely_write_to_database(stmt);

	}
}

} // namesapce
} // namespace
