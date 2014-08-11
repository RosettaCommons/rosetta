// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SecondaryStructureSegmentFeatures.cc
///
/// @brief
/// @author Tim Jacobs

//Unit headers
#include <protocols/features/SecondaryStructureSegmentFeatures.hh>

// Project Headers
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/scoring/dssp/Dssp.hh>

#ifdef WIN32
#include <core/scoring/ScoreFunction.hh>
#endif

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

namespace protocols{
namespace features{

SecondaryStructureSegmentFeatures::SecondaryStructureSegmentFeatures()
{}

SecondaryStructureSegmentFeatures::SecondaryStructureSegmentFeatures(core::scoring::ScoreFunctionOP /*scfxn*/)
{
}

SecondaryStructureSegmentFeatures::SecondaryStructureSegmentFeatures(
		SecondaryStructureSegmentFeatures const & /*src*/ ) : FeaturesReporter(/*src*/)
{
}

SecondaryStructureSegmentFeatures::~SecondaryStructureSegmentFeatures(){}

///@brief return string with class name
std::string
SecondaryStructureSegmentFeatures::type_name() const
{
	return "SecondaryStructureSegmentFeatures";
}

///@brief generate the table schemas and write them to the database
void
SecondaryStructureSegmentFeatures::write_schema_to_db(
	utility::sql_database::sessionOP db_session) const
{
	using namespace basic::database::schema_generator;

	//******secondary structure segments******//
	Column struct_id("struct_id", new DbBigInt(), false);
	Column segment_id("segment_id", new DbInteger(), false);
	Column residue_begin("residue_begin", new DbInteger(), false);
	Column residue_end("residue_end", new DbInteger(), false);
	Column dssp("dssp", new DbText(1), false);

	utility::vector1<Column> pkey_cols;
	pkey_cols.push_back(struct_id);
	pkey_cols.push_back(segment_id);

	Columns foreign_key_columns1;
	foreign_key_columns1.push_back(struct_id);
	utility::vector1< std::string > reference_columns1;
	reference_columns1.push_back("struct_id");
	ForeignKey foreign_key1(foreign_key_columns1, "structures", reference_columns1, true);

	Columns foreign_key_columns2;
	foreign_key_columns2.push_back(struct_id);
	foreign_key_columns2.push_back(residue_begin);
	utility::vector1< std::string > reference_columns2;
	reference_columns2.push_back("struct_id");
	reference_columns2.push_back("resNum");
	ForeignKey foreign_key2(foreign_key_columns2, "residues", reference_columns2, true);

	Columns foreign_key_columns3;
	foreign_key_columns3.push_back(struct_id);
	foreign_key_columns3.push_back(residue_end);
	utility::vector1< std::string > reference_columns3;
	reference_columns3.push_back("struct_id");
	reference_columns3.push_back("resNum");
	ForeignKey foreign_key3(foreign_key_columns3, "residues", reference_columns3, true);

//	Columns foreign_key_columns4;
//	foreign_key_columns4.push_back(dssp);
//	utility::vector1< std::string > reference_columns4;
//	reference_columns4.push_back("code");
//	ForeignKey foreign_key4(foreign_key_columns4, "dssp_codes", reference_columns4, true);

	Schema secondary_structure_segments("secondary_structure_segments", PrimaryKey(pkey_cols));
	secondary_structure_segments.add_foreign_key(foreign_key1);
	secondary_structure_segments.add_foreign_key(foreign_key2);
	secondary_structure_segments.add_foreign_key(foreign_key3);
	secondary_structure_segments.add_column(dssp);
//	secondary_structure_segments.add_foreign_key(foreign_key4);

	secondary_structure_segments.write(db_session);
}

///@brief return the set of features reporters that are required to
///also already be extracted by the time this one is used.
utility::vector1<std::string>
SecondaryStructureSegmentFeatures::features_reporter_dependencies() const
{
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	dependencies.push_back("ResidueSecondaryStructureFeatures");
	return dependencies;
}

void
SecondaryStructureSegmentFeatures::parse_my_tag(
	utility::tag::TagCOP const /*tag*/,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/)
{}

///@brief Return true if the dssp code is irregular (blank), hydrogen
///bonded turn(T) or bend(S)
bool
SecondaryStructureSegmentFeatures::is_loop(std::string dssp_code){
	return
		dssp_code == "T" ||
		dssp_code == "S" ||
		dssp_code == "B" ||
		dssp_code == "I" ||
		dssp_code == " ";
}

bool
SecondaryStructureSegmentFeatures::is_helix(std::string dssp_code){
	return
		dssp_code == "H" ||
		dssp_code == "G";
}

///@brief collect all the feature data for the pose
core::Size
SecondaryStructureSegmentFeatures::report_features(
	core::pose::Pose const & /*pose*/,
	utility::vector1< bool > const & relevant_residues,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session)
{
	using namespace basic::database;

	std::string select_secondary_struct_string=
		"SELECT resNum, dssp\n"
		"FROM residue_secondary_structure\n"
		"WHERE struct_id=?\n"
		"ORDER BY resNum";

	cppdb::statement select_stmt(safely_prepare_statement(select_secondary_struct_string,db_session));
	select_stmt.bind(1,struct_id);
	cppdb::result res=safely_read_from_database(select_stmt);

	std::string sec_structure_statement_string=
	"INSERT INTO secondary_structure_segments (struct_id, segment_id, residue_begin, residue_end,dssp)\n"
	"VALUES (?,?,?,?,?)";

	core::Size segment_counter=0;
	std::string segment_secondary="";
	core::Size segment_begin=1;
	core::Size segment_end;

	// used to determine if segement should be reported based on relevant_residues
	utility::vector1<Size> current_residues;

	core::Size prev_resNum=0;
	core::Size resNum;
	std::string residue_secondary;
	while(res.next()){
	

		res >> resNum >> residue_secondary;

		std::cout << "Current resNum " << resNum << std::endl;

		//Use non-standard 'L' for all loop-like dssp codes
		if(is_loop(residue_secondary))
		{
			residue_secondary="L";
		}
		//include 3/10 helix as 'H', DSSP often classifies ends of helices as 3/10
		else if(is_helix(residue_secondary))
		{
			residue_secondary="H";
		}

		if(residue_secondary != segment_secondary || (resNum != prev_resNum+1 && prev_resNum !=0))
		{
			if(resNum > 1 && check_relevant_residues(relevant_residues, current_residues))
			{
				segment_end=prev_resNum;

				//write segment to DB
				cppdb::statement stmt(basic::database::safely_prepare_statement(sec_structure_statement_string,db_session));
				stmt.bind(1,struct_id);
				stmt.bind(2,segment_counter);
				stmt.bind(3,segment_begin);
				stmt.bind(4,segment_end);
				stmt.bind(5,segment_secondary);
				basic::database::safely_write_to_database(stmt);

				//increment segment secondary structure id
				++segment_counter;
				current_residues.clear();
			}
			segment_secondary = residue_secondary;
			segment_begin=resNum;
		}
		current_residues.push_back(resNum);
		prev_resNum = resNum;
	}

	//write final segment to DB
	cppdb::statement stmt(basic::database::safely_prepare_statement(sec_structure_statement_string,db_session));
	stmt.bind(1,struct_id);
	stmt.bind(2,segment_counter);
	stmt.bind(3,segment_begin);
	stmt.bind(4,resNum);
	stmt.bind(5,segment_secondary);
	basic::database::safely_write_to_database(stmt);

	return 0;
}

} // namespace
} // namespace
