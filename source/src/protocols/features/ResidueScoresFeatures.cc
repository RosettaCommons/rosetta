// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/ResidueScoresFeatures.cc
/// @brief  report residue scores to features Statistics Scientific Benchmark
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/features/ResidueScoresFeatures.hh>
#include <protocols/features/util.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.hh>
#include <core/scoring/ScoreType.hh>
#include <core/types.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>
#include <basic/database/insert_statement_generator/InsertGenerator.hh>
#include <basic/database/insert_statement_generator/RowData.hh>

// Utility Headers
#include <numeric/xyzVector.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <utility/tools/make_vector.hh>
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
#include <protocols/features/ResidueScoresFeaturesCreator.hh>

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
using core::scoring::Energies;
using core::scoring::ScoreFunction;
using core::scoring::ScoreFunctionOP;
using core::scoring::get_score_function;
using core::scoring::ScoreTypeManager;
using core::scoring::ScoreTypes;
using core::scoring::EnergyGraph;
using core::scoring::EnergyEdge;
using core::scoring::LREnergyContainerCOP;
using core::scoring::ResidueNeighborConstIteratorOP;
using core::chemical::aa_vrt;
using protocols::filters::Filters_map;
using basic::datacache::DataMap;
using protocols::moves::Movers_map;
using numeric::xyzVector;
using basic::database::insert_statement_generator::InsertGenerator;
using basic::database::insert_statement_generator::RowData;
using basic::database::insert_statement_generator::RowDataBaseOP;
using utility::tag::TagCOP;
using utility::vector1;
using utility::sql_database::sessionOP;
using utility::tools::make_vector;
using cppdb::statement;

ResidueScoresFeatures::ResidueScoresFeatures() :
	scfxn_(get_score_function())
{}

ResidueScoresFeatures::ResidueScoresFeatures(
	ScoreFunctionOP scfxn) :
	scfxn_(std::move(scfxn))
{}

ResidueScoresFeatures::ResidueScoresFeatures( ResidueScoresFeatures const & src) :
	FeaturesReporter(),
	scfxn_(src.scfxn_)
{}

ResidueScoresFeatures::~ResidueScoresFeatures() = default;

// XRW TEMP string
// XRW TEMP ResidueScoresFeatures::type_name() const { return "ResidueScoresFeatures"; }

void
ResidueScoresFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	write_residue_scores_1b_table_schema(db_session);
	write_residue_scores_2b_table_schema(db_session);
	write_residue_scores_lr_2b_table_schema(db_session);
}

void
ResidueScoresFeatures::write_residue_scores_1b_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column batch_id("batch_id", DbDataTypeOP( new DbInteger() ));
	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column resNum("resNum", DbDataTypeOP( new DbInteger() ));
	Column score_type_id("score_type_id", DbDataTypeOP( new DbInteger() ));
	Column score_value("score_value", DbDataTypeOP( new DbReal() ));
	Column context_dependent("context_dependent", DbDataTypeOP( new DbInteger() ));

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
	write_residue_scores_2b_table_schema_helper( "residue_scores_2b", db_session );
}

void
ResidueScoresFeatures::write_residue_scores_lr_2b_table_schema(
	sessionOP db_session
) const {
	write_residue_scores_2b_table_schema_helper( "residue_scores_lr_2b", db_session );
}

void
ResidueScoresFeatures::write_residue_scores_2b_table_schema_helper(
	std::string name,
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column batch_id("batch_id", DbDataTypeOP( new DbInteger() ));
	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column resNum1("resNum1", DbDataTypeOP( new DbInteger() ));
	Column resNum2("resNum2", DbDataTypeOP( new DbInteger() ));
	Column score_type_id("score_type_id", DbDataTypeOP( new DbInteger() ));
	Column score_value("score_value", DbDataTypeOP( new DbReal() ));
	Column context_dependent("context_dependent", DbDataTypeOP( new DbInteger() ));

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


	Schema table(name, primary_key);
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
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const & /*filters*/,
	Movers_map const & /*movers*/,
	Pose const & /*pose*/
) {
	if ( tag->hasOption("scorefxn") ) {
		string const scorefxn_name(tag->getOption<string>("scorefxn"));
		scfxn_ = data.get_ptr<ScoreFunction>("scorefxns", scorefxn_name);
	}
}

Size
ResidueScoresFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID const struct_id,
	sessionOP db_session
){
	Pose pose_copy(pose);
	(*scfxn_)(pose_copy);
	insert_residue_scores_rows(pose_copy, relevant_residues, struct_id, db_session);

	return 0;
}

/// @details
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
	StructureID const struct_id,
	sessionOP db_session
){

	//calling setup for scoring on a temp copy of the pose, maybe
	//there's a better way of doing this but we can't (and shouldn't)
	//require that the pose be previously setup for scoring before
	//calling this mover
	Pose temp_pose = pose;
	core::pose::symmetry::make_score_function_consistent_with_symmetric_state_of_pose(pose, scfxn_);
	scfxn_->setup_for_scoring(temp_pose);

	Size const batch_id(get_batch_id(struct_id, db_session));


	vector1<bool> relevant_and_virtual_residues(relevant_residues);
	// Since some scores terms, such as elec_dens_fast and constraints,
	// use virtual residues to be compatible with the two-body scoring
	// framework, include virtual residues with the relevant residues so
	// these scores get computed.
	for ( Size i = 1; i <= pose.size(); ++i ) {
		if ( pose.residue( i ).aa() == aa_vrt ) {
			relevant_and_virtual_residues[i] = true;
		}
	}


	insert_one_body_residue_score_rows(
		temp_pose, relevant_and_virtual_residues, batch_id, struct_id, db_session);

	insert_two_body_residue_score_rows(
		temp_pose, relevant_and_virtual_residues, batch_id, struct_id, db_session);

	insert_two_body_long_range_residue_score_rows(
		temp_pose, relevant_and_virtual_residues, batch_id, struct_id, db_session);

} // End function body

void
ResidueScoresFeatures::insert_one_body_residue_score_rows(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size const batch_id,
	StructureID const struct_id,
	sessionOP db_session
) {

	ScoreTypes const & ci_1b( scfxn_->ci_1b_types() );
	ScoreTypes const & cd_1b( scfxn_->cd_1b_types() );

	InsertGenerator insert_onebody("residue_scores_1b");
	insert_onebody.add_column("batch_id");
	insert_onebody.add_column("struct_id");
	insert_onebody.add_column("resNum");
	insert_onebody.add_column("score_type_id");
	insert_onebody.add_column("score_value");
	insert_onebody.add_column("context_dependent");

	RowDataBaseOP batch_id_data( new RowData<Size>("batch_id", batch_id) );
	RowDataBaseOP struct_id_data( new RowData<StructureID>("struct_id", struct_id) );

	EnergyMap emap;

	for ( Size resNum=1; resNum <= pose.size(); ++resNum ) {
		if ( !check_relevant_residues( relevant_residues, resNum ) ) continue;
		Residue const & rsd( pose.residue(resNum) );

		RowDataBaseOP resNum_data( new RowData<Size>("resNum", resNum) );

		{ // Context Independent One Body Energies
			RowDataBaseOP context_dependent_data( new RowData<bool>("context_dependent", false) );

			emap.clear();
			scfxn_->eval_ci_1b(rsd, pose, emap);
			for ( auto st : ci_1b ) {
				if ( !emap[st] ) continue;

				RowDataBaseOP score_type_id_data( new RowData<Size>("score_type_id", st) );
				RowDataBaseOP score_value_data( new RowData<Real>("score_value", emap[st]) );

				insert_onebody.add_row(
					make_vector(
					batch_id_data, struct_id_data, resNum_data,
					score_type_id_data, score_value_data, context_dependent_data));
			}
		}
		{ // Context Dependent One Body Energies
			RowDataBaseOP context_dependent_data( new RowData<bool>("context_dependent", true) );

			emap.clear();
			scfxn_->eval_cd_1b(rsd, pose, emap);
			for ( auto st : cd_1b ) {

				if ( !emap[st] ) continue;

				RowDataBaseOP score_type_id_data( new RowData<Size>("score_type_id", st) );
				RowDataBaseOP score_value_data( new RowData<Real>("score_value", emap[st]) );

				insert_onebody.add_row(
					make_vector(
					batch_id_data, struct_id_data, resNum_data,
					score_type_id_data, score_value_data, context_dependent_data));

			}
		}
	}
	insert_onebody.write_to_database(db_session);
}

void
ResidueScoresFeatures::insert_two_body_residue_score_rows(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size batch_id,
	StructureID const struct_id,
	sessionOP db_session
) {

	// retrieve cached energies object
	Energies const & energies( pose.energies() );
	debug_assert(energies.energies_updated());
	// the neighbor/energy links
	EnergyGraph const & energy_graph( energies.energy_graph() );
	EnergyMap emap;

	ScoreTypes const & ci_2b( scfxn_->ci_2b_types() );
	ScoreTypes const & cd_2b( scfxn_->cd_2b_types() );

	InsertGenerator insert_twobody("residue_scores_2b");
	insert_twobody.add_column("batch_id");
	insert_twobody.add_column("struct_id");
	insert_twobody.add_column("resNum1");
	insert_twobody.add_column("resNum2");
	insert_twobody.add_column("score_type_id");
	insert_twobody.add_column("score_value");
	insert_twobody.add_column("context_dependent");

	RowDataBaseOP batch_id_data( new RowData<Size>("batch_id", batch_id) );
	RowDataBaseOP struct_id_data( new RowData<StructureID>("struct_id", struct_id) );

	for ( Size resNum=1; resNum <= pose.size(); ++resNum ) {
		Residue const & rsd( pose.residue(resNum) );

		// Two Body Energies
		for ( utility::graph::Graph::EdgeListConstIter
				iru  = energy_graph.get_node(resNum)->const_upper_edge_list_begin(),
				irue = energy_graph.get_node(resNum)->const_upper_edge_list_end();
				iru != irue; ++iru ) {
			EnergyEdge const & edge( static_cast< EnergyEdge const &> (**iru) );
			Size const otherResNum( edge.get_second_node_ind() );

			if ( !check_relevant_residues( relevant_residues, resNum, otherResNum ) ) continue;

			Size resNum1, resNum2;
			if ( resNum < otherResNum ) {
				resNum1 = resNum;
				resNum2 = otherResNum;
			} else {
				resNum1 = otherResNum;
				resNum2 = resNum;
			}

			Residue const & otherRsd( pose.residue(otherResNum) );

			RowDataBaseOP resNum1_data( new RowData<Size>("resNum1", resNum1) );
			RowDataBaseOP resNum2_data( new RowData<Size>("resNum2", resNum2) );

			{ // Context Independent Two Body Energies

				RowDataBaseOP context_dependent_data( new RowData<bool>("context_dependent", false) );

				emap.clear();
				scfxn_->eval_ci_2b(rsd, otherRsd, pose, emap);
				for ( auto st : ci_2b ) {
					if ( !emap[st] ) continue;

					RowDataBaseOP score_type_id_data( new RowData<Size>("score_type_id", st) );
					RowDataBaseOP score_value_data( new RowData<Real>("score_value", emap[st]) );

					insert_twobody.add_row(
						make_vector(
						batch_id_data, struct_id_data, resNum1_data, resNum2_data,
						score_type_id_data, score_value_data, context_dependent_data));

				}
			}
			{ // Context Dependent Two Body Energies
				RowDataBaseOP context_dependent_data( new RowData<bool>("context_dependent", true) );

				EnergyMap emap;
				scfxn_->eval_cd_2b(rsd, otherRsd, pose, emap);
				for ( auto st : cd_2b ) {
					if ( !emap[st] ) continue;

					RowDataBaseOP score_type_id_data( new RowData<Size>("score_type_id", st) );
					RowDataBaseOP score_value_data( new RowData<Real>("score_value", emap[st]) );

					insert_twobody.add_row(
						make_vector(
						batch_id_data, struct_id_data, resNum1_data, resNum2_data,
						score_type_id_data, score_value_data, context_dependent_data));
				}
			}
		}
	}
	insert_twobody.write_to_database(db_session);
}


void
ResidueScoresFeatures::insert_two_body_long_range_residue_score_rows(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size batch_id,
	StructureID const struct_id,
	sessionOP db_session
) {


	ScoreTypes const & ci_lr_2b( scfxn_->ci_lr_2b_types() );
	ScoreTypes const & cd_lr_2b( scfxn_->cd_lr_2b_types() );

	InsertGenerator insert_twobody_longrange("residue_scores_lr_2b");
	insert_twobody_longrange.add_column("batch_id");
	insert_twobody_longrange.add_column("struct_id");
	insert_twobody_longrange.add_column("resNum1");
	insert_twobody_longrange.add_column("resNum2");
	insert_twobody_longrange.add_column("score_type_id");
	insert_twobody_longrange.add_column("score_value");
	insert_twobody_longrange.add_column("context_dependent");

	RowDataBaseOP batch_id_data( new RowData<Size>("batch_id", batch_id) );
	RowDataBaseOP struct_id_data( new RowData<StructureID>("struct_id", struct_id) );

	EnergyMap emap;

	{ // Context Independent Long Range Two Body Energies
		RowDataBaseOP context_dependent_data( new RowData<bool>("context_dependent", false) );

		for ( auto
				iter = scfxn_->ci_lr_2b_methods_begin(),
				iter_end = scfxn_->ci_lr_2b_methods_end();
				iter != iter_end; ++iter ) {
			LREnergyContainerCOP lrec =
				pose.energies().long_range_container((*iter)->long_range_type());
			if ( !lrec || lrec->empty() ) continue; // only score non-emtpy energies.

			// Potentially O(N^2) operation...
			for ( Size resNum = 1; resNum <= pose.size(); ++resNum ) {

				for ( ResidueNeighborConstIteratorOP
						rni = lrec->const_upper_neighbor_iterator_begin( resNum ),
						rniend = lrec->const_upper_neighbor_iterator_end( resNum );
						(*rni) != (*rniend); ++(*rni) ) {
					Size const otherResNum(rni->upper_neighbor_id());
					if ( !check_relevant_residues( relevant_residues, resNum, otherResNum ) ) continue;

					Size resNum1, resNum2;
					if ( resNum < otherResNum ) {
						resNum1 = resNum;
						resNum2 = otherResNum;
					} else {
						resNum1 = otherResNum;
						resNum2 = resNum;
					}

					debug_assert(rni->energy_computed());
					emap.zero();
					rni->retrieve_energy( emap );

					RowDataBaseOP resNum1_data( new RowData<Size>("resNum1", resNum1) );
					RowDataBaseOP resNum2_data( new RowData<Size>("resNum2", resNum2) );

					for ( auto st : ci_lr_2b ) {
						if ( !emap[st] ) continue;

						RowDataBaseOP score_type_id_data( new RowData<Size>("score_type_id", st) );
						RowDataBaseOP score_value_data( new RowData<Real>("score_value", emap[st]) );

						insert_twobody_longrange.add_row(
							make_vector(
							batch_id_data, struct_id_data, resNum1_data, resNum2_data,
							score_type_id_data, score_value_data, context_dependent_data));
					}
				}
			}
		}
	}

	/////////////////////////////////////////////////////
	///  Context Dependent Long Range Two Body methods
	{
		RowDataBaseOP context_dependent_data( new RowData<bool>("context_dependent", true) );

		for ( auto
				iter = scfxn_->cd_lr_2b_methods_begin(),
				iter_end = scfxn_->cd_lr_2b_methods_end();
				iter != iter_end; ++iter ) {
			LREnergyContainerCOP lrec =
				pose.energies().long_range_container((*iter)->long_range_type());
			if ( !lrec || lrec->empty() ) continue; // only score non-emtpy energies.

			// Potentially O(N^2) operation...
			for ( Size resNum = 1; resNum <= pose.size(); ++resNum ) {
				if ( !relevant_residues[resNum] ) continue;

				for ( ResidueNeighborConstIteratorOP
						rni = lrec->const_upper_neighbor_iterator_begin( resNum ),
						rniend = lrec->const_upper_neighbor_iterator_end( resNum );
						(*rni) != (*rniend); ++(*rni) ) {
					Size const otherResNum(rni->upper_neighbor_id());
					if ( !relevant_residues[otherResNum] ) continue;

					Size resNum1, resNum2;
					if ( resNum < otherResNum ) {
						resNum1 = resNum;
						resNum2 = otherResNum;
					} else {
						resNum1 = otherResNum;
						resNum2 = resNum;
					}

					debug_assert(rni->energy_computed());
					emap.zero();
					rni->retrieve_energy( emap );

					RowDataBaseOP resNum1_data( new RowData<Size>("resNum1", resNum1) );
					RowDataBaseOP resNum2_data( new RowData<Size>("resNum2", resNum2) );

					for (
							auto
							st = ci_lr_2b.begin(), ste = cd_lr_2b.end();
							st != ste; ++st ) {
						if ( !emap[*st] ) continue;

						RowDataBaseOP score_type_id_data( new RowData<Size>("score_type_id", *st) );
						RowDataBaseOP score_value_data( new RowData<Real>("score_value", emap[*st]) );

						insert_twobody_longrange.add_row(
							make_vector(
							batch_id_data, struct_id_data, resNum1_data, resNum2_data,
							score_type_id_data, score_value_data, context_dependent_data));
					}
				}
			}
		}
	}
	insert_twobody_longrange.write_to_database(db_session);
}

std::string ResidueScoresFeatures::type_name() const {
	return class_name();
}

std::string ResidueScoresFeatures::class_name() {
	return "ResidueScoresFeatures";
}

void ResidueScoresFeatures::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "scorefxn", xs_string, "Scorefunction to be used" );

	protocols::features::xsd_type_definition_w_attributes( xsd, class_name(), "Record all per-residue and per-residue-pair scores given a scoring function", attlist );
}

std::string ResidueScoresFeaturesCreator::type_name() const {
	return ResidueScoresFeatures::class_name();
}

protocols::features::FeaturesReporterOP
ResidueScoresFeaturesCreator::create_features_reporter() const {
	return protocols::features::FeaturesReporterOP( new ResidueScoresFeatures );
}

void ResidueScoresFeaturesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ResidueScoresFeatures::provide_xml_schema( xsd );
}



} // namesapce
} // namespace
