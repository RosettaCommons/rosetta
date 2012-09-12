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
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/features/ResidueScoresFeatures.hh>
#include <protocols/features/util.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/types.hh>
#include <protocols/moves/DataMap.hh>
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>

// Utility Headers
#include <numeric/xyzVector.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

// External Headers
#include <cppdb/frontend.h>
#include <boost/uuid/uuid_io.hpp>

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

void
ResidueScoresFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	write_residue_scores_1b_table_schema(db_session);
	write_residue_scores_2b_table_schema(db_session);
}

void
ResidueScoresFeatures::write_residue_scores_1b_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column batch_id("batch_id", new DbInteger());
	Column struct_id("struct_id", new DbUUID());
	Column resNum("resNum", new DbInteger());
	Column score_type_id("score_type_id", new DbInteger());
	Column score_value("score_value", new DbReal());
	Column context_dependent("context_dependent", new DbInteger());

	Columns primary_key_columns;
	primary_key_columns.push_back(batch_id);
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(resNum);
	primary_key_columns.push_back(score_type_id);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns1;
	foreign_key_columns1.push_back(struct_id);
	foreign_key_columns1.push_back(resNum);
	vector1< std::string > reference_columns1;
	reference_columns1.push_back("struct_id");
	reference_columns1.push_back("resNum");
	ForeignKey foreign_key1(foreign_key_columns1, "residues", reference_columns1, true);

	Columns foreign_key_columns2;
	foreign_key_columns2.push_back(batch_id);
	foreign_key_columns2.push_back(score_type_id);
	vector1< std::string > reference_columns2;
	reference_columns2.push_back("batch_id");
	reference_columns2.push_back("score_type_id");
	ForeignKey foreign_key2(foreign_key_columns2, "score_types", reference_columns2, true);


	Schema table("residue_scores_1b", primary_key);
	table.add_foreign_key(foreign_key1);
	table.add_foreign_key(foreign_key2);
	table.add_column(score_value);
	table.add_column(context_dependent);

	table.write(db_session);
}

void
ResidueScoresFeatures::write_residue_scores_2b_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column batch_id("batch_id", new DbInteger());
	Column struct_id("struct_id", new DbUUID());
	Column resNum1("resNum1", new DbInteger());
	Column resNum2("resNum2", new DbInteger());
	Column score_type_id("score_type_id", new DbInteger());
	Column score_value("score_value", new DbReal());
	Column context_dependent("context_dependent", new DbInteger());

	Columns primary_key_columns;
	primary_key_columns.push_back(batch_id);
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(resNum1);
	primary_key_columns.push_back(resNum2);
	primary_key_columns.push_back(score_type_id);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns1;
	foreign_key_columns1.push_back(struct_id);
	foreign_key_columns1.push_back(resNum1);
	vector1< std::string > reference_columns1;
	reference_columns1.push_back("struct_id");
	reference_columns1.push_back("resNum");
	ForeignKey foreign_key1(foreign_key_columns1, "residues", reference_columns1, true);

	Columns foreign_key_columns2;
	foreign_key_columns2.push_back(struct_id);
	foreign_key_columns2.push_back(resNum2);
	vector1< std::string > reference_columns2;
	reference_columns2.push_back("struct_id");
	reference_columns2.push_back("resNum");
	ForeignKey foreign_key2(foreign_key_columns2, "residues", reference_columns2, true);

	Columns foreign_key_columns3;
	foreign_key_columns3.push_back(batch_id);
	foreign_key_columns3.push_back(score_type_id);
	vector1< std::string > reference_columns3;
	reference_columns3.push_back("batch_id");
	reference_columns3.push_back("score_type_id");
	ForeignKey foreign_key3(foreign_key_columns3, "score_types", reference_columns3, true);


	Schema table("residue_scores_2b", primary_key);
	table.add_foreign_key(foreign_key1);
	table.add_foreign_key(foreign_key2);
	table.add_foreign_key(foreign_key3);
	table.add_column(score_value);
	table.add_column(context_dependent);

	table.write(db_session);
}

utility::vector1<std::string>
ResidueScoresFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	dependencies.push_back("ScoreTypeFeatures");
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
	boost::uuids::uuid const struct_id,
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
	boost::uuids::uuid const struct_id,
	sessionOP db_session
){

	// I would like to assert that this has been called, but I don't know how
	//scfxn_.setup_for_scoring( pose );

	//calling setup for scoring on a temp copy of the pose, maybe there's a better way of doing this
	//but we can't (and shouldn't) require that the pose be previously setup for scoring before calling this mover
	Pose temp_pose = pose;
	scfxn_->setup_for_scoring(temp_pose);

	Size const batch_id(get_batch_id(struct_id, db_session));

	ScoreTypes ci_1b( scfxn_->ci_1b_types() );
	ScoreTypes cd_1b( scfxn_->cd_1b_types() );
	ScoreTypes ci_2b( scfxn_->ci_2b_types() );
	ScoreTypes cd_2b( scfxn_->cd_2b_types() );

	std::string oneb_string = "INSERT INTO residue_scores_1b (batch_id, struct_id, resNum, score_type_id, score_value, context_dependent) VALUES (?,?,?,?,?,?);";
	std::string twob_string = "INSERT INTO residue_scores_2b (batch_id, struct_id, resNum1, resNum2, score_type_id, score_value, context_dependent) VALUES (?,?,?,?,?,?,?);";

	statement oneb_stmt(basic::database::safely_prepare_statement(oneb_string,db_session));
	statement twob_stmt(basic::database::safely_prepare_statement(twob_string,db_session));

	for(Size resNum=1; resNum <= temp_pose.total_residue(); ++resNum){
		if(!relevant_residues[resNum]) continue;
		Residue rsd( temp_pose.residue(resNum) );
		{ // Context Independent One Body Energies
			EnergyMap emap;
			scfxn_->eval_ci_1b(rsd, temp_pose, emap);
			for(ScoreTypes::const_iterator st = ci_1b.begin(), ste = ci_1b.end(); st != ste; ++st){
				if(!emap[*st]) continue;

				Real const score_value( emap[*st] );
				bool const context_dependent(false);

				oneb_stmt.bind(1, batch_id);
				oneb_stmt.bind(2, struct_id);
				oneb_stmt.bind(3, resNum);
				oneb_stmt.bind(4, *st);
				oneb_stmt.bind(5, score_value);
				oneb_stmt.bind(6, context_dependent);
				basic::database::safely_write_to_database(oneb_stmt);
			}
		}
		{ // Context Dependent One Body Energies
			EnergyMap emap;
			scfxn_->eval_cd_1b(rsd, temp_pose, emap);
			for(ScoreTypes::const_iterator
				st = cd_1b.begin(), ste = cd_1b.end();
				st != ste; ++st){

				if(!emap[*st]) continue;

				Real const score_value( emap[*st] );
				bool const context_dependent(true);

				oneb_stmt.bind(1, batch_id);
				oneb_stmt.bind(2, struct_id);
				oneb_stmt.bind(3, resNum);
				oneb_stmt.bind(4, *st);
				oneb_stmt.bind(5, score_value);
				oneb_stmt.bind(6, context_dependent);
				basic::database::safely_write_to_database(oneb_stmt);

			}
		}

		// Two Body Energies
		for(Size otherResNum=resNum+1; otherResNum <= temp_pose.total_residue(); ++otherResNum){
			if(!relevant_residues[otherResNum]) continue;
			if(!scfxn_->are_they_neighbors(temp_pose, resNum, otherResNum)) continue;
			Residue otherRsd( temp_pose.residue(otherResNum) );
			{ // Context Independent Two Body Energies
				EnergyMap emap;
				scfxn_->eval_ci_2b(rsd, otherRsd, temp_pose, emap);
				for(ScoreTypes::const_iterator st = ci_2b.begin(), ste = ci_2b.end(); st != ste; ++st){
					if(!emap[*st]) continue;

					Real const score_value( emap[*st] );
					bool const context_dependent(false);

					twob_stmt.bind(1, batch_id);
					twob_stmt.bind(2, struct_id);
					twob_stmt.bind(3, resNum);
					twob_stmt.bind(4, otherResNum);
					twob_stmt.bind(5, *st);
					twob_stmt.bind(6, score_value);
					twob_stmt.bind(7, context_dependent);
					basic::database::safely_write_to_database(twob_stmt);

				}
			}
			{ // Context Dependent Two Body Energies
				EnergyMap emap;
				scfxn_->eval_cd_2b(rsd, otherRsd, temp_pose, emap);
				for(ScoreTypes::const_iterator st = cd_2b.begin(), ste = cd_2b.end(); st != ste; ++st){
					if(!emap[*st]) continue;

					Real const score_value( emap[*st] );
					bool const context_dependent(true);

					twob_stmt.bind(1, batch_id);
					twob_stmt.bind(2, struct_id);
					twob_stmt.bind(3, resNum);
					twob_stmt.bind(4, otherResNum);
					twob_stmt.bind(5, *st);
					twob_stmt.bind(6, score_value);
					twob_stmt.bind(7, context_dependent);
					basic::database::safely_write_to_database(twob_stmt);

				}
			}
		} // End Two Body Energies
	} // End res1 for loop
} // End function body


} // namesapce
} // namespace
