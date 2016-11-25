// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/ResidueTotalScoresFeatures.cc
/// @brief  report the total per residue score to a features database
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/features/ResidueTotalScoresFeatures.hh>
#include <protocols/features/util.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/types.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/database/sql_utils.hh>
#include <basic/Tracer.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>

// Utility Headers
#include <numeric/xyzVector.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <cmath>
#include <utility/excn/Exceptions.hh>
#include <sstream>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/features/feature_schemas.hh>
#include <protocols/features/ResidueTotalScoresFeaturesCreator.hh>

namespace protocols {
namespace features {

using std::endl;
using std::string;
using std::stringstream;
using core::Size;
using core::Real;
using core::pose::Pose;
using core::pose::PoseOP;
using core::conformation::Residue;
using core::scoring::EnergyMap;
using core::scoring::ScoreFunction;
using core::scoring::ScoreFunctionOP;
using core::scoring::get_score_function;
using core::scoring::ScoreTypeManager;
using core::scoring::ScoreTypes;
using core::pose::symmetry::make_score_function_consistent_with_symmetric_state_of_pose;
using protocols::filters::Filters_map;
using basic::datacache::DataMap;
using protocols::moves::Movers_map;
using numeric::xyzVector;
using utility::tag::TagCOP;
using utility::vector1;
using basic::database::safely_write_to_database;
using utility::sql_database::sessionOP;
using cppdb::statement;

static THREAD_LOCAL basic::Tracer TR( "protocols.features.ResidueTotalScoresFeatures" );

ResidueTotalScoresFeatures::ResidueTotalScoresFeatures() :
	scfxn_(get_score_function())
{}

ResidueTotalScoresFeatures::ResidueTotalScoresFeatures(
	ScoreFunctionOP scfxn) :
	scfxn_(std::move(scfxn))
{}

ResidueTotalScoresFeatures::ResidueTotalScoresFeatures( ResidueTotalScoresFeatures const & src) :
	FeaturesReporter(),
	scfxn_(src.scfxn_)
{}

ResidueTotalScoresFeatures::~ResidueTotalScoresFeatures()
= default;

// XRW TEMP string
// XRW TEMP ResidueTotalScoresFeatures::type_name() const { return "ResidueTotalScoresFeatures"; }

void
ResidueTotalScoresFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	write_residue_total_scores_table_schema(db_session);
}

void
ResidueTotalScoresFeatures::write_residue_total_scores_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column resNum("resNum", DbDataTypeOP( new DbInteger() ));
	Column score_value("score_value", DbDataTypeOP( new DbReal() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(resNum);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(struct_id);
	foreign_key_columns.push_back(resNum);
	vector1< string > reference_columns;
	reference_columns.push_back("struct_id");
	reference_columns.push_back("resNum");
	ForeignKey foreign_key(foreign_key_columns, "residues", reference_columns, true);

	Schema table("residue_total_scores", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(score_value);

	table.write(db_session);
}

utility::vector1<string>
ResidueTotalScoresFeatures::features_reporter_dependencies() const {
	utility::vector1<string> dependencies;
	dependencies.push_back("ResidueFeatures");
	dependencies.push_back("ScoreTypeFeatures");
	return dependencies;
}

void
ResidueTotalScoresFeatures::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const & /*filters*/,
	Movers_map const & /*movers*/,
	Pose const & /*pose*/
) {
	if ( tag->hasOption("scorefxn") ) {
		string const scorefxn_name(tag->getOption<string>("scorefxn"));
		scfxn_ = data.get_ptr<ScoreFunction>("scorefxns", scorefxn_name);
	} else {
		stringstream error_msg;
		error_msg
			<< "The " << type_name() << " reporter requires a 'scorefxn' tag:" << endl
			<< endl
			<< "    <feature name=" << type_name() <<" scorefxn=(name_of_score_function) />" << endl;
		throw utility::excn::EXCN_RosettaScriptsOption(error_msg.str());
	}
}

Size
ResidueTotalScoresFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID const struct_id,
	sessionOP db_session
){
	core::pose::symmetry::make_score_function_consistent_with_symmetric_state_of_pose(pose, scfxn_);
	insert_residue_total_scores_rows(pose, relevant_residues, struct_id, db_session );

	return 0;
}

/// @details
void
ResidueTotalScoresFeatures::insert_residue_total_scores_rows(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID const struct_id,
	sessionOP db_session
) const {
	if ( !scfxn_->energy_method_options().hbond_options().decompose_bb_hb_into_pair_energies() ) {
		TR.Warning << "The backbone-backbone hydrogen bonds being stored in the per-residue scores. Please enable the decompose_bb_hb_into_pair_energies option if you want this behavior." << endl;
	}


	//calling setup for scoring on a temp copy of the pose, maybe
	//there's a better way of doing this but we can't (and shouldn't)
	//require that the pose be scored before calling this.
	Pose temp_pose = pose;
	(*scfxn_)(temp_pose);

	//Size const batch_id(get_batch_id(struct_id, db_session));

	string const stmt_string("INSERT INTO residue_total_scores (struct_id, resNum, score_value ) VALUES (?,?,?);");

	statement stmt(
		basic::database::safely_prepare_statement(stmt_string, db_session));

	for ( Size resNum=1; resNum <= temp_pose.size(); ++resNum ) {
		if ( !check_relevant_residues( relevant_residues, resNum ) ) continue;

		Real const score_value(temp_pose.energies().residue_total_energy(resNum));

		stmt.bind(1, struct_id);
		stmt.bind(2, resNum);
		stmt.bind(3, score_value);
		safely_write_to_database(stmt);

	} // End rsd for loop


} // End function body

std::string ResidueTotalScoresFeatures::type_name() const {
	return class_name();
}

std::string ResidueTotalScoresFeatures::class_name() {
	return "ResidueTotalScoresFeatures";
}

void ResidueTotalScoresFeatures::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute(
		"scorefxn", xs_string,
		"score function name");

	protocols::features::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"report the total per residue score to a features database",
		attlist );
}

std::string ResidueTotalScoresFeaturesCreator::type_name() const {
	return ResidueTotalScoresFeatures::class_name();
}

protocols::features::FeaturesReporterOP
ResidueTotalScoresFeaturesCreator::create_features_reporter() const {
	return protocols::features::FeaturesReporterOP( new ResidueTotalScoresFeatures );
}

void ResidueTotalScoresFeaturesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ResidueTotalScoresFeatures::provide_xml_schema( xsd );
}



} // namesapce
} // namespace
