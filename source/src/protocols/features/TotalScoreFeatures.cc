// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Unit headers
#include <protocols/features/TotalScoreFeatures.hh>
#include <protocols/features/TotalScoreFeaturesCreator.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Utility headers
#include <basic/database/insert_statement_generator/InsertGenerator.hh>
#include <basic/database/insert_statement_generator/RowData.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/sql_utils.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/tools/make_vector.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

// RosettaScripts headers
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.hh>

// C++ headers
#include <string>
#include <sstream>

namespace protocols {
namespace features {

using namespace std;
using core::Size;
using core::Real;
using core::pose::Pose;
using core::scoring::ScoreFunction;
using core::scoring::ScoreFunctionOP;
using core::scoring::ScoreFunctionCOP;
using utility::tools::make_vector;

FeaturesReporterOP TotalScoreFeaturesCreator::create_features_reporter() const {
	return FeaturesReporterOP( new TotalScoreFeatures );
}

std::string TotalScoreFeaturesCreator::type_name() const {
	return "TotalScoreFeatures";
}

TotalScoreFeatures::TotalScoreFeatures()
: scorefxn_(core::scoring::get_score_function()) {}

TotalScoreFeatures::TotalScoreFeatures(ScoreFunctionOP scorefxn)
: scorefxn_(std::move(scorefxn)) {}

TotalScoreFeatures::~TotalScoreFeatures() = default;

string TotalScoreFeatures::type_name() const {
	return "TotalScoreFeatures";
}

ScoreFunctionCOP TotalScoreFeatures::scorefxn() const {
	return scorefxn_;
}

void TotalScoreFeatures::scorefxn(ScoreFunctionOP scorefxn) {
	scorefxn_ = scorefxn;
}

void TotalScoreFeatures::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/) {

	if ( tag->hasOption("scorefxn") ) {
		string scorefxn_name = tag->getOption<string>("scorefxn");
		scorefxn_ = data.get_ptr<ScoreFunction>("scorefxns", scorefxn_name);
	} else {
		stringstream error_msg;
		error_msg
			<< "The " << type_name() << " reporter requires a 'scorefxn' tag:" << endl
			<< endl
			<< "    <feature name=" << type_name() <<" scorefxn=(name_of_score_function) />" << endl;
		throw utility::excn::EXCN_RosettaScriptsOption(error_msg.str());
	}
}

void TotalScoreFeatures::write_schema_to_db(
	utility::sql_database::sessionOP db_session) const {

	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt ), false);
	Column score("score", DbDataTypeOP( new DbReal ), false);

	Schema total_scores("total_scores", PrimaryKey(struct_id));

	total_scores.add_column(score);
	total_scores.write(db_session);
}

Size TotalScoreFeatures::report_features(
	core::pose::Pose const & pose,
	utility::vector1<bool> const & /*relevant_residues*/,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session) {

	// Calculate the total score of the pose.  The entire pose is included in the
	// score calculation; the relevant residues_parameter is ignored.

	Pose non_const_pose = pose;
	runtime_assert(scorefxn_.get() != nullptr);
	core::pose::symmetry::make_score_function_consistent_with_symmetric_state_of_pose(pose, scorefxn_);
	Real total_score = scorefxn_->score(non_const_pose);

	// Write the total score to the database.

	using namespace basic::database::insert_statement_generator;

	InsertGenerator insert_statement("total_scores");
	insert_statement.add_column("struct_id");
	insert_statement.add_column("score");
	insert_statement.add_row(make_vector<RowDataBaseOP>(
		RowDataBaseOP( new RowData<StructureID>("struct_id", struct_id) ),
		RowDataBaseOP( new RowData<Real>("score", total_score) )));

	insert_statement.write_to_database(db_session);

	return 0;
}


}
}
