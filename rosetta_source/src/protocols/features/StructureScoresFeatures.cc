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
#include <core/graph/Graph.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <sstream>
#include <string>

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
using core::scoring::hbonds::HBondSet;
using core::scoring::hbonds::get_hbond_energies;
using core::scoring::EnergiesCacheableDataType::HBOND_SET;
using core::Size;
using core::Real;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;

StructureScoresFeatures::StructureScoresFeatures() :
	initialized_score_types_(false),
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
	initialized_score_types_(src.initialized_score_types_),
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
			"CREATE TABLE IF NOT EXISTS score_types (\n"
			"	protocol_id INTEGER,\n"
			"	score_type_id INTEGER PRIMARY KEY,\n"
			"	score_type_name TEXT,\n"
			"	FOREIGN KEY (protocol_id)\n"
			"		REFERENCES protocols (protocol_id)\n"
			"		DEFERRABLE INITIALLY DEFERRED);\n"
			"\n"
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
			"CREATE TABLE IF NOT EXISTS score_types (\n"
			"	protocol_id INTEGER,\n"
			"	score_type_id INTEGER PRIMARY KEY,\n"
			"	score_type_name TEXT,\n"
			"	FOREIGN KEY (protocol_id) REFERENCES protocols (protocol_id));\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS structure_scores (\n"
			"	struct_id INTEGER,\n"
			"	score_type_id INTEGER,\n"
			"	score_value INTEGER,\n"
			"	FOREIGN KEY (struct_id) REFERENCES structures (struct_id),\n"
			"	FOREIGN KEY (score_type_id) REFERENCES score_types (score_type_id),\n"
			"	PRIMARY KEY (struct_id, score_type_id));";
	}else
	{
		return "";
	}

}

Size
StructureScoresFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size struct_id,
	sessionOP db_session
){
	insert_score_type_rows(struct_id, db_session);
	insert_structure_score_rows(pose, relevant_residues, struct_id, db_session);

	return 0;
}

void
StructureScoresFeatures::insert_score_type_rows(
	Size struct_id,
	sessionOP db_session
) {
	std::string db_mode(basic::options::option[basic::options::OptionKeys::inout::database_mode]);
	Size protocol_id;
	result res = (*db_session) <<
		"SELECT protocol_id FROM structures WHERE struct_id=?" << struct_id;
	if(res.next()){
		res >> protocol_id;
	} else {
		stringstream error_message;
		error_message << "Structure with struct_id '" << struct_id << "' does not exist in the database." << endl;
		utility_exit_with_message(error_message.str());
	}

	//cppdb::transaction transact_guard(*db_session);
	for(Size score_type_id=1; score_type_id <= n_score_types; ++score_type_id){
		ScoreType type(static_cast<ScoreType>(score_type_id));

		string const score_type( ScoreTypeManager::name_from_score_type(type) );

		// I <3 all the tiny differences between SQL dialects,
		//it makes the code so much more vibrant and diverse
		if(db_mode == "sqlite3")
		{
			statement stmt = (*db_session)
				<< "INSERT OR IGNORE INTO score_types VALUES (?,?,?);"
				<< protocol_id
				<< score_type_id
				<< score_type;
			stmt.exec();
		}else if(db_mode == "mysql")
		{
			statement stmt = (*db_session)
				<< "INSERT IGNORE INTO score_types VALUES (?,?,?);"
				<< protocol_id
				<< score_type_id
				<< score_type;
			stmt.exec();
		}else
		{
			utility_exit_with_message("the database mode needs to be 'mysql' or 'sqlite3'");
		}

	}
	//transact_guard.commit();
	initialized_score_types_ = true;
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


	for(Size score_type_id=1; score_type_id <= n_score_types; ++score_type_id){
		ScoreType type(static_cast<ScoreType>(score_type_id));
		Real const score_value( energies.weights()[type] * emap[type] );
		if(!score_value) continue;
		statement stmt = (*db_session)
			<< "INSERT INTO structure_scores VALUES (?,?,?);"
			<< struct_id
			<< score_type_id
			<< score_value;
		stmt.exec();
	}
}

} // namesapce
} // namespace
