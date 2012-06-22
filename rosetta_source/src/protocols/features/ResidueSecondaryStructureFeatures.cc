// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ResidueSecondaryStructureFeatures.cc
/// @brief  report ResidueSecondaryStructure geometry and scores to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/ResidueSecondaryStructureFeatures.hh>

//External
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/assign/list_of.hpp>

// Project Headers
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>
#include <core/scoring/dssp/Dssp.hh>

//Basic Headers
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/Constraint.hh>

// Utility Headers
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>

// External Headers
#include <cppdb/frontend.h>

namespace protocols{
namespace features{

using std::string;
using core::scoring::dssp::Dssp;
using core::Size;
using core::conformation::Residue;
using core::pose::Pose;
using core::scoring::ScoreFunctionOP;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using basic::Tracer;

static Tracer TR("protocols.features.ResidueSecondaryStructureFeatures");

ResidueSecondaryStructureFeatures::ResidueSecondaryStructureFeatures() {}

ResidueSecondaryStructureFeatures::ResidueSecondaryStructureFeatures(ResidueSecondaryStructureFeatures const &) :
	FeaturesReporter()
{}

ResidueSecondaryStructureFeatures::~ResidueSecondaryStructureFeatures() {}

string
ResidueSecondaryStructureFeatures::type_name() const { return "ResidueSecondaryStructureFeatures"; }

void
ResidueSecondaryStructureFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const{
	using namespace basic::database::schema_generator;
	using namespace basic::database;
	using namespace boost::assign;

	Column code("code", DbText(1), false);
	Column label("label", DbText(), false);
	Schema dssp_codes("dssp_codes", PrimaryKey(code));
	dssp_codes.add_column(label);

	dssp_codes.write(db_session);

	std::vector<std::string> dssp_cols = list_of("code")("label");
	insert_or_ignore("dssp_codes", dssp_cols, list_of("'H'")("'H: a-Helix'"), db_session);
	insert_or_ignore("dssp_codes", dssp_cols, list_of("'E'")("'E: b-Sheet'"), db_session);
	insert_or_ignore("dssp_codes", dssp_cols, list_of("'T'")("'T: HB Turn'"), db_session);
	insert_or_ignore("dssp_codes", dssp_cols, list_of("'G'")("'G: 3/10 Helix'"), db_session);
	insert_or_ignore("dssp_codes", dssp_cols, list_of("'B'")("'B: b-Bridge'"), db_session);
	insert_or_ignore("dssp_codes", dssp_cols, list_of("'S'")("'S: Bend'"), db_session);
	insert_or_ignore("dssp_codes", dssp_cols, list_of("'I'")("'I: pi-Helix'"), db_session);
	insert_or_ignore("dssp_codes", dssp_cols, list_of("' '")("'Irregular'"), db_session);

	/******residue_secondary_structure******/
	Column struct_id("struct_id",DbUUID(), false);
	Column resNum("resNum",DbInteger(), false);
	Column dssp("dssp",DbText(1));

	utility::vector1<Column> sec_struct_pkey_cols;
	sec_struct_pkey_cols.push_back(struct_id);
	sec_struct_pkey_cols.push_back(resNum);

	utility::vector1<Column> fkey_cols;
	fkey_cols.push_back(struct_id);
	fkey_cols.push_back(resNum);

	utility::vector1<std::string> fkey_reference_cols;
	fkey_reference_cols.push_back("struct_id");
	fkey_reference_cols.push_back("resNum");

	ForeignKey dssp_fk(dssp, "dssp_codes", "code", true /*defer*/);

	Schema residue_secondary_structure("residue_secondary_structure", PrimaryKey(sec_struct_pkey_cols));

	residue_secondary_structure.add_column(struct_id);
	residue_secondary_structure.add_column(resNum);
	residue_secondary_structure.add_column(dssp);

	residue_secondary_structure.add_foreign_key(ForeignKey(fkey_cols, "residues", fkey_reference_cols, true));
	residue_secondary_structure.add_foreign_key(dssp_fk);

	residue_secondary_structure.write(db_session);

	/******helix_segments******/
	Column helix_id("helix_id",DbInteger(), false);
	Column residue_begin("residue_begin",DbInteger(), false);
	Column residue_end("residue_end",DbInteger(), false);

	utility::vector1<Column> helix_pkey_cols;
	helix_pkey_cols.push_back(struct_id);
	helix_pkey_cols.push_back(helix_id);

	utility::vector1<Column> fkey_cols_begin;
	fkey_cols_begin.push_back(struct_id);
	fkey_cols_begin.push_back(residue_begin);

	utility::vector1<Column> fkey_cols_end;
	fkey_cols_end.push_back(struct_id);
	fkey_cols_end.push_back(residue_end);

	Schema helix_segments("helix_segments", PrimaryKey(helix_pkey_cols));

	helix_segments.add_column(struct_id);
	helix_segments.add_column(helix_id);
	helix_segments.add_column(residue_begin);
	helix_segments.add_column(residue_end);

	helix_segments.add_foreign_key(ForeignKey(fkey_cols_begin, "residues", fkey_reference_cols, true));
	helix_segments.add_foreign_key(ForeignKey(fkey_cols_end, "residues", fkey_reference_cols, true));

	helix_segments.write(db_session);

	/******beta_segments******/
	Column beta_id("beta_id",DbInteger(), false);

	utility::vector1<Column> beta_pkey_cols;
	beta_pkey_cols.push_back(struct_id);
	beta_pkey_cols.push_back(beta_id);

	Schema beta_segments("beta_segments", PrimaryKey(beta_pkey_cols));

	beta_segments.add_column(struct_id);
	beta_segments.add_column(beta_id);
	beta_segments.add_column(residue_begin);
	beta_segments.add_column(residue_end);

	beta_segments.add_foreign_key(ForeignKey(fkey_cols_begin, "residues", fkey_reference_cols, true));
	beta_segments.add_foreign_key(ForeignKey(fkey_cols_end, "residues", fkey_reference_cols, true));

	beta_segments.write(db_session);
}

utility::vector1<std::string>
ResidueSecondaryStructureFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	return dependencies;
}

Size
ResidueSecondaryStructureFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	boost::uuids::uuid struct_id,
	sessionOP db_session
){
	// compute dssp
	core::scoring::dssp::Dssp all_dssp(pose);

	//stores the secondary structure for the current stretch of secondary structure
	string segment_secondary;
	Size segment_begin=1;
	Size segment_end;
	Size helix_counter=1;
	Size beta_counter=1;
	//Size loop_counter=1;

	//Create the statement strings outside the loops so we don't need to rcreate them for every residue
	std::string sec_structure_statement_string = "INSERT INTO residue_secondary_structure (struct_id, resNum, dssp) VALUES (?,?,?);";
	std::string helix_segment_statement_string = "INSERT INTO helix_segments (struct_id, helix_id, residue_begin, residue_end) VALUES (?,?,?,?);";
	std::string beta_segment_statement_string = "INSERT INTO beta_segments (struct_id, beta_id, residue_begin, residue_end) VALUES (?,?,?,?);";
	for(Size resNum=1; resNum <= pose.total_residue(); ++resNum){
		if(!relevant_residues[resNum]) continue;


		if(!pose.residue(resNum).is_protein()){
			// Due to limitations with the current DSSP code,
			// the indexing gets off after a non-protein residue
			// and leads to a segmentation fault.
			break;
		}

		string residue_secondary = string(1, all_dssp.get_dssp_secstruct(resNum));
		if(residue_secondary != segment_secondary){

			if(resNum > 1){
				segment_end=resNum-1;
				if(segment_secondary == "H"){

					statement stmt(basic::database::safely_prepare_statement(helix_segment_statement_string,db_session));
					stmt.bind(1,struct_id);
					stmt.bind(2,helix_counter);
					stmt.bind(3,segment_begin);
					stmt.bind(4,segment_end);

					basic::database::safely_write_to_database(stmt);

					++helix_counter;
				}
				else if(segment_secondary == "E"){

					statement stmt(basic::database::safely_prepare_statement(beta_segment_statement_string,db_session));
					stmt.bind(1,struct_id);
					stmt.bind(2,beta_counter);
					stmt.bind(3,segment_begin);
					stmt.bind(4,segment_end);

					basic::database::safely_write_to_database(stmt);

					++beta_counter;
				}
			}
			segment_secondary = residue_secondary;
			segment_begin=resNum;
		}

		statement stmt(basic::database::safely_prepare_statement(sec_structure_statement_string,db_session));
		stmt.bind(1,struct_id);
		stmt.bind(2,resNum);
		stmt.bind(3,residue_secondary);
		basic::database::safely_write_to_database(stmt);

	}
	return 0;
}


} // namesapce
} // namespace
