// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/StructureScoresFeatures.cc
/// @brief  report protocol level features to features statistics scientific benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/StructureScoresFeatures.hh>

// Platform Headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <protocols/moves/DataMap.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <basic/database/sql_utils.hh>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <sstream>
#include <string>

#include <core/scoring/ScoreFunction.hh>
#include <utility/vector0.hh>


namespace protocols{
namespace features{

using std::string;
using std::stringstream;
using std::endl;
using core::pose::Pose;
using core::graph::Graph;
using core::scoring::Energies;
using core::scoring::EnergyGraph;
using core::scoring::EnergyMap;
using core::scoring::EnergyEdge;
using core::scoring::getScoreFunction;
using core::scoring::hbond_sr_bb;
using core::scoring::hbond_lr_bb;
using core::scoring::n_score_types;
using core::scoring::ScoreTypeManager;
using core::scoring::ScoreType;
using core::scoring::ScoreFunctionOP;
using core::scoring::ScoreFunction;
using core::scoring::hbonds::HBondSet;
using core::scoring::hbonds::get_hbond_energies;
using core::scoring::EnergiesCacheableDataType::HBOND_SET;
using core::Size;
using core::Real;
using protocols::filters::Filters_map;
using protocols::moves::DataMap;
using protocols::moves::Movers_map;
using utility::vector1;
using utility::sql_database::sessionOP;
using utility::tag::TagPtr;
using cppdb::statement;
using cppdb::result;

StructureScoresFeatures::StructureScoresFeatures() :
	scfxn_(getScoreFunction())
{}

StructureScoresFeatures::StructureScoresFeatures(
	ScoreFunctionOP scfxn
) :
	scfxn_(scfxn)
{
	if ( scfxn_ == 0 ) {
		utility_exit_with_message( "StructureScoresFeatures may not be constructed with a null-pointer ScoreFunctionOP" );
	}
}

StructureScoresFeatures::StructureScoresFeatures(
	StructureScoresFeatures const & src
) :
	FeaturesReporter(),
	scfxn_(src.scfxn_)
{}

StructureScoresFeatures::~StructureScoresFeatures() {}

string
StructureScoresFeatures::type_name() const { return "StructureScoresFeatures"; }

string
StructureScoresFeatures::schema() const {
	std::string db_mode(basic::options::option[basic::options::OptionKeys::inout::database_mode]);

	if(db_mode == "sqlite3")
	{
		return
			"CREATE TABLE IF NOT EXISTS structure_scores (\n"
			"	struct_id INTEGER,\n"
			"	score_type_id INTEGER,\n"
			"	score_value INTEGER,\n"
			"	FOREIGN KEY (struct_id)\n"
			"		REFERENCES structures (struct_id)\n"
			"		DEFERRABLE INITIALLY DEFERRED,\n"
			"	FOREIGN KEY (score_type_id)\n"
			"		REFERENCES score_types (score_type_id)\n"
			"		DEFERRABLE INITIALLY DEFERRED,\n"
			"	PRIMARY KEY (struct_id, score_type_id));";
	}else if(db_mode == "mysql")
	{
		return
			"CREATE TABLE IF NOT EXISTS structure_scores (\n"
			"	struct_id BIGINT UNSIGNED REFERENCES structures (struct_id),\n"
			"	score_type_id INTEGER REFERENCES score_types (score_type_id),\n"
			"	score_value INTEGER,\n"
			"	PRIMARY KEY (struct_id, score_type_id));\n";
	}else
	{
		return "";
	}

}

utility::vector1<std::string>
StructureScoresFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("StructureFeatures");
	dependencies.push_back("ScoreTypeFeatures");
	return dependencies;
}

void
StructureScoresFeatures::parse_my_tag(
	TagPtr const tag,
	DataMap & data,
	Filters_map const & /*filters*/,
	Movers_map const & /*movers*/,
	Pose const & /*pose*/
) {
	if(tag->hasOption("scorefxn")){
		string scorefxn_name = tag->getOption<string>("scorefxn");
		scfxn_ = data.get<ScoreFunction*>("scorefxns", scorefxn_name);
	} else {
		stringstream error_msg;
		error_msg
			<< "The " << type_name() << " reporter requires a 'scorefxn' tag:" << endl
			<< endl
			<< "    <feature name=" << type_name() <<" scorefxn=(name_of_score_function) />" << endl;
		utility_exit_with_message(error_msg.str());
	}
}





Size
StructureScoresFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size struct_id,
	sessionOP db_session
){
	if(!pose.energies().energies_updated()){
		stringstream err_msg;
		err_msg << "Attempting to extract structure score features from pose, however the energies are not up to date. Please score the pose with a ScoreFunction first.";
		utility_exit_with_message(err_msg.str());
	}

	insert_structure_score_rows(pose, relevant_residues, struct_id, db_session);
	return 0;
}

void StructureScoresFeatures::delete_record(
	Size struct_id,
	utility::sql_database::sessionOP db_session
){

	std::string statement_string = "DELETE FROM structure_scores WHERE struct_id = ?;\n";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);
	basic::database::safely_write_to_database(stmt);

}

void
StructureScoresFeatures::insert_structure_score_rows(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size struct_id,
	sessionOP db_session
) const {

	Energies const & energies(pose.energies());
	EnergyMap emap;
	scfxn_->get_sub_score(pose, relevant_residues, emap);

	core::Real total_score= 0.0;

	std::string statement_string = "INSERT INTO structure_scores VALUES (?,?,?);";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	for(Size score_type_id=1; score_type_id <= n_score_types; ++score_type_id){
		ScoreType type(static_cast<ScoreType>(score_type_id));
		Real const score_value( energies.weights()[type] * emap[type] );
		if(!score_value) continue;
		total_score += score_value;
		stmt.bind(1,struct_id);
		stmt.bind(2,score_type_id);
		stmt.bind(3,score_value);
		basic::database::safely_write_to_database(stmt);
	}
	// add the total_score type
	stmt.bind(1,struct_id);
	stmt.bind(2,n_score_types);
	stmt.bind(3,total_score);
	basic::database::safely_write_to_database(stmt);
}

} // namesapce
} // namespace
