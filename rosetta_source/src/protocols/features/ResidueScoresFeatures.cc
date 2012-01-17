// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ResidueScoresFeatures.cc
/// @brief  report residue scores to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/ResidueScoresFeatures.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/types.hh>
#include <protocols/moves/DataMap.hh>
#include <basic/database/sql_utils.hh>

// Utility Headers
#include <numeric/xyzVector.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <cmath>
#include <sstream>

namespace protocols{
namespace features{

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
using core::scoring::getScoreFunction;
using core::scoring::ScoreTypeManager;
using core::scoring::ScoreTypes;
using protocols::filters::Filters_map;
using protocols::moves::DataMap;
using protocols::moves::Movers_map;
using numeric::xyzVector;
using utility::tag::TagPtr;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;


ResidueScoresFeatures::ResidueScoresFeatures() :
	scfxn_(getScoreFunction())
{}

ResidueScoresFeatures::ResidueScoresFeatures(
	ScoreFunctionOP scfxn) :
	scfxn_(scfxn)
{}

ResidueScoresFeatures::ResidueScoresFeatures( ResidueScoresFeatures const & src) :
	FeaturesReporter(),
	scfxn_(src.scfxn_)
{}

ResidueScoresFeatures::~ResidueScoresFeatures()
{}

string
ResidueScoresFeatures::type_name() const { return "ResidueScoresFeatures"; }

string
ResidueScoresFeatures::schema() const {
	return
		"CREATE TABLE IF NOT EXISTS residue_scores_1b (\n"
		"	struct_id INTEGER,\n"
		"	resNum INTEGER,\n"
		"	score_type TEXT,\n"
		"	score_value REAL,\n"
		"	context_dependent  INTEGER,\n"
		"	FOREIGN KEY (struct_id, resNum)\n"
		"		REFERENCES residues (struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	PRIMARY KEY(struct_id, resNum, score_type));\n"
		"\n"
		"CREATE TABLE IF NOT EXISTS residue_scores_2b (\n"
		"	struct_id INTEGER,\n"
		"	resNum1 INTEGER,\n"
		"	resNum2 INTEGER,\n"
		"	score_type TEXT,\n"
		"	score_value REAL,\n"
		"	context_dependent  INTEGER,\n"
		"	FOREIGN KEY (struct_id, resNum1)\n"
		"		REFERENCES residues (struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	FOREIGN KEY (struct_id, resNum2)\n"
		"		REFERENCES residues (struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	PRIMARY KEY(struct_id, resNum1, resNum2, score_type));\n";
}

utility::vector1<std::string>
ResidueScoresFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	return dependencies;
}

void
ResidueScoresFeatures::parse_my_tag(
	TagPtr const tag,
	DataMap & data,
	Filters_map const & /*filters*/,
	Movers_map const & /*movers*/,
	Pose const & /*pose*/
) {
	if(tag->hasOption("scorefxn")){
		string const scorefxn_name(tag->getOption<string>("scorefxn"));
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
ResidueScoresFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size const struct_id,
	sessionOP db_session
){
	insert_residue_scores_rows(pose, relevant_residues, struct_id, db_session );

	return 0;
}

///@details
///
/// * Score types are either one body, two body, or whole structure and
/// can either dependend on the context or not.
///
/// * Whole structure terms are reported in the structure_scores table
/// along with the totals from the remaining terms.
///
/// * The one and two body terms are broken up into two different
/// tables because they are parametrized differently (one residue vs
/// two residues).
///
/// * Residues are identified by Rosetta's residue numbering scheme.
/// To convert to the PDB residue numbering scheme, join with the
/// residues_pdb table.
///
/// * Although two body terms can be with in the same residue, these
/// 'intrares' score terms are reported with the two body terms where
/// resNum == otherResNum.
///
/// * Two body terms always are reported so that resNum1 <= resNum2.
///
/// * Values for score terms are only reported if they are non-zero.
///
/// * Values of two body energies are only reported when both residues
/// are true in the relevant_residues vector.

void
ResidueScoresFeatures::insert_residue_scores_rows(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size const struct_id,
	sessionOP db_session
){

	// I would like to assert that this has been called, but I don't know how
	//scfxn_.setup_for_scoring( pose );

	ScoreTypes ci_1b( scfxn_->ci_1b_types() );
	ScoreTypes cd_1b( scfxn_->cd_1b_types() );
	ScoreTypes ci_2b( scfxn_->ci_2b_types() );
	ScoreTypes cd_2b( scfxn_->cd_2b_types() );

	std::string oneb_string = "INSERT INTO residue_scores_1b VALUES (?,?,?,?,?);";
	std::string twob_string = "INSERT INTO residue_scores_2b VALUES (?,?,?,?,?,?);";

	statement oneb_stmt(basic::database::safely_prepare_statement(oneb_string,db_session));
	statement twob_stmt(basic::database::safely_prepare_statement(twob_string,db_session));

	for(Size resNum=1; resNum <= pose.total_residue(); ++resNum){
		if(!relevant_residues[resNum]) continue;
		Residue rsd( pose.residue(resNum) );
		{ // Context Independent One Body Energies
			EnergyMap emap;
			scfxn_->eval_ci_1b(rsd, pose, emap);
			for(ScoreTypes::const_iterator st = ci_1b.begin(), ste = ci_1b.end(); st != ste; ++st){
				if(!emap[*st]) continue;

				string const score_type( ScoreTypeManager::name_from_score_type(*st) );
				Real const score_value( emap[*st] );

				oneb_stmt.bind(1,struct_id);
				oneb_stmt.bind(2,resNum);
				oneb_stmt.bind(3,score_type);
				oneb_stmt.bind(4,score_value);
				oneb_stmt.bind(5,false);
				basic::database::safely_write_to_database(oneb_stmt);
			}
		}
		{ // Context Dependent One Body Energies
			EnergyMap emap;
			scfxn_->eval_cd_1b(rsd, pose, emap);
			for(ScoreTypes::const_iterator
				st = cd_1b.begin(), ste = cd_1b.end();
				st != ste; ++st){

				if(!emap[*st]) continue;

				string const score_type( ScoreTypeManager::name_from_score_type(*st) );
				Real const score_value( emap[*st] );

				oneb_stmt.bind(1,struct_id);
				oneb_stmt.bind(2,resNum);
				oneb_stmt.bind(3,score_type);
				oneb_stmt.bind(4,score_value);
				oneb_stmt.bind(5,true);
				basic::database::safely_write_to_database(oneb_stmt);

			}
		}

		// Two Body Energies
		for(Size otherResNum=resNum+1; otherResNum <= pose.total_residue(); ++otherResNum){
			if(!relevant_residues[otherResNum]) continue;
			if(!scfxn_->are_they_neighbors(pose, resNum, otherResNum)) continue;
			Residue otherRsd( pose.residue(otherResNum) );
			{ // Context Independent Two Body Energies
				EnergyMap emap;
				scfxn_->eval_ci_2b(rsd, otherRsd, pose, emap);
				for(ScoreTypes::const_iterator st = ci_2b.begin(), ste = ci_2b.end(); st != ste; ++st){
					if(!emap[*st]) continue;

					string const score_type(ScoreTypeManager::name_from_score_type(*st));
					Real const score_value( emap[*st] );

					twob_stmt.bind(1,struct_id);
					twob_stmt.bind(2,resNum);
					twob_stmt.bind(3,otherResNum);
					twob_stmt.bind(4,score_type);
					twob_stmt.bind(5,score_value);
					twob_stmt.bind(6,false);
					basic::database::safely_write_to_database(twob_stmt);

				}
			}
			{ // Context Dependent Two Body Energies
				EnergyMap emap;
				scfxn_->eval_cd_2b(rsd, otherRsd, pose, emap);
				for(ScoreTypes::const_iterator st = cd_2b.begin(), ste = cd_2b.end(); st != ste; ++st){
					if(!emap[*st]) continue;
					string const score_type(ScoreTypeManager::name_from_score_type(*st));
					Real const score_value( emap[*st] );

					twob_stmt.bind(1,struct_id);
					twob_stmt.bind(2,resNum);
					twob_stmt.bind(3,otherResNum);
					twob_stmt.bind(4,score_type);
					twob_stmt.bind(5,score_value);
					twob_stmt.bind(6,true);
					basic::database::safely_write_to_database(twob_stmt);

				}
			}
		} // End Two Body Energies
	} // End res1 for loop
} // End function body


} // namesapce
} // namespace
