// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/StructureScoresFeatures.cc
/// @brief  report structure score features to a features database
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/StructureScoresFeatures.hh>
#include <protocols/features/util.hh>

//External
#include <cppdb/frontend.h>

//Basic Headers
#include <basic/Tracer.hh>
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/Constraint.hh>

// Platform Headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <basic/database/sql_utils.hh>
#include <utility/tools/make_vector.hh>

#include <basic/database/insert_statement_generator/InsertGenerator.hh>
#include <basic/database/insert_statement_generator/RowData.hh>

// C++ Headers
#include <sstream>
#include <string>

#include <core/scoring/ScoreFunction.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector0.hh>


namespace protocols{
namespace features{

static thread_local basic::Tracer TR( "protocols.features.StructureScoresFeatures" );

using std::string;
using std::stringstream;
using std::endl;
using core::pose::Pose;
using core::graph::Graph;
using core::scoring::Energies;
using core::scoring::EnergyGraph;
using core::scoring::EnergyMap;
using core::scoring::EnergyEdge;
using core::scoring::get_score_function;
using core::scoring::hbond_sr_bb;
using core::scoring::hbond_lr_bb;
using core::scoring::n_score_types;
using core::scoring::total_score;
using core::scoring::ScoreTypeManager;
using core::scoring::ScoreType;
using core::scoring::ScoreFunctionOP;
using core::scoring::ScoreFunction;
using core::chemical::aa_vrt;
using core::scoring::hbonds::HBondSet;
using core::scoring::hbonds::get_hbond_energies;
using core::scoring::EnergiesCacheableDataType::HBOND_SET;
using core::Size;
using core::Real;
using protocols::filters::Filters_map;
using basic::datacache::DataMap;
using protocols::moves::Movers_map;
using utility::vector1;
using utility::sql_database::sessionOP;
using utility::tag::TagCOP;
using cppdb::statement;
using cppdb::result;
using basic::database::insert_statement_generator::InsertGenerator;
using basic::database::insert_statement_generator::RowDataBaseOP;
using basic::database::insert_statement_generator::RowData;

StructureScoresFeatures::StructureScoresFeatures() :
	scfxn_(get_score_function())
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

void
StructureScoresFeatures::write_schema_to_db(
	utility::sql_database::sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	//******structure_scores******//
	Column batch_id("batch_id", DbDataTypeOP( new DbInteger() ), false);
	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ), false);
	Column score_type_id("score_type_id", DbDataTypeOP( new DbInteger() ), false);
	Column score_value("score_value", DbDataTypeOP( new DbReal() ), false);

	utility::vector1<Column> pkey_cols;
	pkey_cols.push_back(batch_id);
	pkey_cols.push_back(struct_id);
	pkey_cols.push_back(score_type_id);

	Columns foreign_key_columns1;
	foreign_key_columns1.push_back(struct_id);
	vector1< std::string > reference_columns1;
	reference_columns1.push_back("struct_id");
	ForeignKey foreign_key1(foreign_key_columns1, "structures", reference_columns1, true);


	Columns foreign_key_columns2;
	foreign_key_columns2.push_back(batch_id);
	foreign_key_columns2.push_back(score_type_id);
	vector1< std::string > reference_columns2;
	reference_columns2.push_back("batch_id");
	reference_columns2.push_back("score_type_id");
	ForeignKey foreign_key2(foreign_key_columns2, "score_types", reference_columns2, true);

	Schema structure_scores("structure_scores", PrimaryKey(pkey_cols));
	structure_scores.add_column(score_value);

	structure_scores.add_foreign_key(foreign_key1);
	structure_scores.add_foreign_key(foreign_key2);

	structure_scores.write(db_session);

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
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const & /*filters*/,
	Movers_map const & /*movers*/,
	Pose const & /*pose*/
) {
	if(tag->hasOption("scorefxn")){
		string scorefxn_name = tag->getOption<string>("scorefxn");
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
StructureScoresFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID struct_id,
	sessionOP db_session
){
	core::pose::symmetry::make_score_function_consistent_with_symmetric_state_of_pose(pose, scfxn_);
	insert_structure_score_rows(pose, relevant_residues, struct_id, db_session);
	return 0;
}

void StructureScoresFeatures::delete_record(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session
){

	std::string statement_string = "DELETE FROM structure_scores WHERE struct_id = ?;\n";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);
	basic::database::safely_write_to_database(stmt);

}

void
StructureScoresFeatures::compute_energies(
	Pose const & pose_in,
	vector1<bool> const & relevant_residues,
	EnergyMap & emap
) const {
	// We only get const access to the pose, so make a copy of it so it
	// can be scored. TODO: This is a bit wasteful, as we don't actually
	// need a copy of the coodinates part of the pose only the
	// energies part of the pose.
	Pose pose = pose_in;
	(*scfxn_)(pose);

	vector1<bool> relevant_and_virtual_residues(relevant_residues);

	// if relevant residues includes all of the residues, then include
	// the whole structure scores, otherwise just the one-body and
	// two-body scores.
	bool all_residues(true);

	if(relevant_residues_mode_ == RelevantResiduesMode::Inclusive){
		TR.Warning << "StructureScoresFeatures is currently not compatible with inclusive RelevantResiduesMode::Inclusive, treating as Explicit." << std::endl;
	}

	// Since some scores terms, such as elec_dens_fast and constraints,
	// use virtual residues to be compatible with the two-body scoring
	// framework, include virtual residues with the relevant residues so
	// these scores get computed.
	for(Size i = 1; i <= pose.total_residue(); ++i){
		if(pose.residue( i ).aa() == aa_vrt){
			relevant_and_virtual_residues[i] = true;
		} else if (relevant_and_virtual_residues[i] == false){
			all_residues = false;
		}
	}
	if(all_residues){
		emap = pose.energies().total_energies();
	} else {
		scfxn_->get_sub_score(pose, relevant_and_virtual_residues, emap);
	}
}


void
StructureScoresFeatures::insert_structure_score_rows(
	Pose const & pose_in,
	vector1< bool > const & relevant_residues,
	StructureID struct_id,
	sessionOP db_session
) const {

	InsertGenerator structure_scores_insert("structure_scores");
	structure_scores_insert.add_column("batch_id");
	structure_scores_insert.add_column("struct_id");
	structure_scores_insert.add_column("score_type_id");
	structure_scores_insert.add_column("score_value");

	EnergyMap emap;
	compute_energies(pose_in, relevant_residues, emap);

	core::Real total_score_value= 0.0;
	Size const batch_id(get_batch_id(struct_id, db_session));
	RowDataBaseOP struct_id_data( new RowData<StructureID>("struct_id",struct_id) );
	RowDataBaseOP batch_id_data( new RowData<Size>("batch_id",batch_id) );

	for(Size score_type_id=1; score_type_id <= n_score_types; ++score_type_id){
		ScoreType type(static_cast<ScoreType>(score_type_id));
		Real const score_value( (*scfxn_)[type] * emap[type] );
		if(!score_value) continue;
		total_score_value += score_value;

		RowDataBaseOP score_type_id_data( new RowData<Size>("score_type_id",score_type_id) );
		RowDataBaseOP score_value_data( new RowData<Real>("score_value",score_value) );

		structure_scores_insert.add_row(
			utility::tools::make_vector(batch_id_data,struct_id_data,score_type_id_data,score_value_data));

	}

	// add the total_score type
	RowDataBaseOP total_id_data( new RowData<Size>("score_type_id",total_score) );
	RowDataBaseOP total_score_data( new RowData<Real>("score_value",total_score_value) );

	structure_scores_insert.add_row(
		utility::tools::make_vector(batch_id_data,struct_id_data,total_id_data,total_score_data));
	structure_scores_insert.write_to_database(db_session);
	}

} // namesapce
} // namespace
