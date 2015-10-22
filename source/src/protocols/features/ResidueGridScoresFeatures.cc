// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ResidueGridScoresFeatures.cc
/// @brief detailed per atom scores of Scoring Grids
/// @author Sam DeLuca

#include <protocols/features/ResidueGridScoresFeatures.hh>
#include <protocols/features/ResidueGridScoresFeaturesCreator.hh>
#include <protocols/qsar/scoring_grid/GridManager.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>
#include <basic/database/sql_utils.hh>

#include <basic/database/insert_statement_generator/InsertGenerator.hh>
#include <basic/database/insert_statement_generator/RowData.hh>

#include <utility/tools/make_vector.hh>
#include <utility/tag/Tag.hh>

#include <utility/excn/Exceptions.hh>
#include <core/pose/util.hh>

namespace protocols {
namespace features {

ResidueGridScoresFeaturesCreator::ResidueGridScoresFeaturesCreator()
{

}

ResidueGridScoresFeaturesCreator::~ResidueGridScoresFeaturesCreator()
{

}

protocols::features::FeaturesReporterOP ResidueGridScoresFeaturesCreator::create_features_reporter() const
{
	return protocols::features::FeaturesReporterOP( new ResidueGridScoresFeatures );
}

std::string ResidueGridScoresFeaturesCreator::type_name() const
{
	return "ResidueGridScoresFeatures";
}

ResidueGridScoresFeatures::ResidueGridScoresFeatures() : chain_(' ')
{

}

ResidueGridScoresFeatures::ResidueGridScoresFeatures(ResidueGridScoresFeatures const & src) : FeaturesReporter(src),
	chain_(src.chain_)
{

}


ResidueGridScoresFeatures::~ResidueGridScoresFeatures()
{

}

std::string ResidueGridScoresFeatures::type_name() const
{
	return "ResidueGridScoresFeatures";
}

void ResidueGridScoresFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const
{
	using namespace basic::database::schema_generator;
	Column struct_id("struct_id",DbDataTypeOP( new DbBigInt() ));
	Column grid_name("grid_name",DbDataTypeOP( new DbTextKey() ));
	Column seqpos("seqpos",DbDataTypeOP( new DbInteger() ));
	Column atomno("atomno",DbDataTypeOP( new DbInteger() ));
	Column score("score",DbDataTypeOP( new DbReal() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(grid_name);
	primary_key_columns.push_back(seqpos);
	primary_key_columns.push_back(atomno);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(struct_id);
	foreign_key_columns.push_back(seqpos);
	foreign_key_columns.push_back(atomno);

	utility::vector1< std::string > reference_columns;
	reference_columns.push_back("struct_id");
	reference_columns.push_back("seqpos");
	reference_columns.push_back("atomno");
	ForeignKey foreign_key(foreign_key_columns,"residue_atom_coords",reference_columns,true);

	Schema table("residue_grid_scores",primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(score);

	table.write(db_session);


}

utility::vector1<std::string> ResidueGridScoresFeatures::features_reporter_dependencies() const
{
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueConformationFeatures");
	return dependencies;
}

core::Size ResidueGridScoresFeatures::report_features(
	core::pose::Pose const & pose,
	utility::vector1< bool > const & relevant_residues,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session)
{
	using basic::database::insert_statement_generator::InsertGenerator;
	using basic::database::insert_statement_generator::RowDataBaseOP;
	using basic::database::insert_statement_generator::RowData;

	protocols::qsar::scoring_grid::GridManager* grid_manager = qsar::scoring_grid::GridManager::get_instance();

	if ( grid_manager->size()==0 ) {
		utility_exit_with_message("In order to use the ResidueGridScoresFeatures reporter you must define at least one scoring grid");
	}

	core::Size chain_id = core::pose::get_chain_id_from_chain(chain_,pose);
	core::Size jump_id = core::pose::get_jump_id_from_chain(chain_,pose);
	core::Vector const center(protocols::geometry::downstream_centroid_by_jump(pose,jump_id));
	grid_manager->initialize_all_grids(center);
	grid_manager->update_grids(pose,center);


	InsertGenerator grid_insert("residue_grid_scores");
	grid_insert.add_column("struct_id");
	grid_insert.add_column("grid_name");
	grid_insert.add_column("seqpos");
	grid_insert.add_column("atomno");
	grid_insert.add_column("score");

	RowDataBaseOP struct_id_data( new RowData<StructureID>("struct_id",struct_id) );
	for ( Size i = pose.conformation().chain_begin(chain_id); i <= pose.conformation().chain_end(chain_id); ++i ) {

		if ( !check_relevant_residues(relevant_residues, i) ) continue;

		RowDataBaseOP seqpos_data( new RowData<core::Size>("seqpos",i) );

		core::conformation::Residue const & residue(pose.residue(i));
		for ( Size atomno = 1; atomno <= residue.natoms(); ++atomno ) {
			RowDataBaseOP atomno_data( new RowData<core::Size>("atomno",atomno) );
			std::map<std::string,core::Real> atom_map = grid_manager->atom_score(pose,residue,atomno);
			for ( std::map<std::string,core::Real>::const_iterator score_it = atom_map.begin(); score_it != atom_map.end(); ++score_it ) {
				RowDataBaseOP grid_name_data( new RowData<std::string>("grid_name",score_it->first) );
				RowDataBaseOP score_data( new RowData<core::Real>("score",score_it->second) );
				grid_insert.add_row(
					utility::tools::make_vector(struct_id_data,seqpos_data,atomno_data,grid_name_data,score_data));
			}
		}
	}
	grid_insert.write_to_database(db_session);
	return 0;
}

void ResidueGridScoresFeatures::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & ,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/)
{
	if ( !tag->hasOption("chain") ) {
		throw utility::excn::EXCN_RosettaScriptsOption("The ResidueGridScoresFeatures reporter requires a Chain tag");

	}

	chain_ = tag->getOption<char>("chain");

}

}
}
