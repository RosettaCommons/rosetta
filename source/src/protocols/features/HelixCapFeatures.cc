// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file HelixCapFeatures.cc
///
/// @brief
/// @author Ben Borgo
/// @author Jim Havranek

//Unit headers
#include <protocols/features/HelixCapFeatures.hh>

// Project Headers
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <string>
#include <cppdb/frontend.h>


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

namespace protocols {
namespace features {

HelixCapFeatures::HelixCapFeatures()
{}

HelixCapFeatures::HelixCapFeatures(core::scoring::ScoreFunctionOP /*scfxn*/ )
{
}

HelixCapFeatures::HelixCapFeatures( HelixCapFeatures const & /*src*/ ) : FeaturesReporter(/*src*/)
{
}

HelixCapFeatures::~HelixCapFeatures(){}

/// @brief return string with class name
std::string
HelixCapFeatures::type_name() const
{
	return "HelixCapFeatures";
}

/// @brief generate the table schemas and write them to the database
void
HelixCapFeatures::write_schema_to_db(
	utility::sql_database::sessionOP db_session) const
{
	using namespace basic::database::schema_generator;

	//******secondary structure segments******//
	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ), false);
	Column resnum("resnum", DbDataTypeOP( new DbInteger() ), false);
	Column type("type", DbDataTypeOP( new DbText(1) ), false);
	Column im4("im4", DbDataTypeOP( new DbInteger() ), false);
	Column im3("im3", DbDataTypeOP( new DbInteger() ), false);
	Column im2("im2", DbDataTypeOP( new DbInteger() ), false);
	Column im1("im1", DbDataTypeOP( new DbInteger() ), false);
	Column ip1("ip1", DbDataTypeOP( new DbInteger() ), false);
	Column ip2("ip2", DbDataTypeOP( new DbInteger() ), false);
	Column ip3("ip3", DbDataTypeOP( new DbInteger() ), false);
	Column ip4("ip4", DbDataTypeOP( new DbInteger() ), false);
	Column dssp("dssp", DbDataTypeOP( new DbText(1) ), false);

	Columns pkey_cols;
	pkey_cols.push_back(struct_id);
	pkey_cols.push_back(resnum);
	PrimaryKey primary_key( pkey_cols );

	Columns foreign_key_columns1;
	foreign_key_columns1.push_back(struct_id);
	utility::vector1< std::string > reference_columns1;
	reference_columns1.push_back("struct_id");
	ForeignKey foreign_key1(foreign_key_columns1, "structures", reference_columns1, true);

	Columns foreign_key_columns2;
	foreign_key_columns2.push_back(struct_id);
	foreign_key_columns2.push_back(resnum);
	utility::vector1< std::string > reference_columns2;
	reference_columns2.push_back("struct_id");
	reference_columns2.push_back("resnum");
	ForeignKey foreign_key2(foreign_key_columns2, "residues", reference_columns2, true);

	Schema secondary_structure_punctuation("secondary_structure_punctuation", primary_key);
	secondary_structure_punctuation.add_foreign_key(foreign_key1);
	secondary_structure_punctuation.add_foreign_key(foreign_key2);
	secondary_structure_punctuation.add_column(dssp);
	secondary_structure_punctuation.add_column(type);
	secondary_structure_punctuation.add_column(im4);
	secondary_structure_punctuation.add_column(im3);
	secondary_structure_punctuation.add_column(im2);
	secondary_structure_punctuation.add_column(im1);
	secondary_structure_punctuation.add_column(ip1);
	secondary_structure_punctuation.add_column(ip2);
	secondary_structure_punctuation.add_column(ip3);
	secondary_structure_punctuation.add_column(ip4);

	secondary_structure_punctuation.write(db_session);
}

/// @brief return the set of features reporters that are required to
///also already be extracted by the time this one is used.
utility::vector1<std::string>
HelixCapFeatures::features_reporter_dependencies() const
{
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	dependencies.push_back("ResidueSecondaryStructureFeatures");
	return dependencies;
}

void
HelixCapFeatures::parse_my_tag(
	utility::tag::TagCOP /*tag*/,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/)
{}

/// @brief Return true if the dssp code is irregular (blank), hydrogen
///bonded turn(T) or bend(S)
bool
HelixCapFeatures::is_loop(std::string dssp_code){
	return
		dssp_code == "T" ||
		dssp_code == "S" ||
		dssp_code == "B" ||
		dssp_code == "G" ||
		dssp_code == "I" ||
		dssp_code == " ";
}

/// @brief collect all the feature data for the pose
core::Size
HelixCapFeatures::report_features(
	core::pose::Pose const & pose,
	utility::vector1< bool > const & /*relevant_residues*/,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session)
{
	using namespace basic::database;

	std::string sec_structure_statement_string=
		"INSERT INTO secondary_structure_punctuation (struct_id, resnum, dssp, type, im4, im3, im2, im1, ip1, ip2, ip3, ip4 ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)";

	cppdb::statement stmt(basic::database::safely_prepare_statement(sec_structure_statement_string,db_session));

	core::scoring::dssp::Dssp all_dssp(pose);

	for ( Size i = 5; i <= pose.total_residue()-5; ++i ) {

		core::conformation::Residue const & resi = pose.residue(i);
		if ( !resi.is_protein() ) continue;

		std::string ss (1, all_dssp.get_dssp_secstruct(i));
		std::string ssC (1, all_dssp.get_dssp_secstruct(i+1));
		std::string ssN (1, all_dssp.get_dssp_secstruct(i-1));

		if ( (ss == "G") || (ss == "I") ) {
			ss = "H";
		}
		if ( (ssC == "G") || (ssC == "I") ) {
			ssC = "H";
		}
		if ( (ssN == "G") || (ssN == "I") ) {
			ssN = "H";
		}

		if ( (ss != ssC) && ((ss == "H") || (ss == "E")) ) {
			std::string ct = "cterm";

			stmt.bind(1,struct_id);
			stmt.bind(2,i);
			stmt.bind(3,ss);
			stmt.bind(4,ct);
			stmt.bind(5,i-4);
			stmt.bind(6,i-3);
			stmt.bind(7,i-2);
			stmt.bind(8,i-1);
			stmt.bind(9,i+1);
			stmt.bind(10,i+2);
			stmt.bind(11,i+3);
			stmt.bind(12,i+4);
			basic::database::safely_write_to_database(stmt);
		} else if ( (ss != ssN) && ((ss == "H") || (ss == "E")) ) {
			// Get the helix length
			Size work_resid( i+1 );
			Size helix_length( 1 );
			while ( work_resid <= pose.total_residue() && all_dssp.get_dssp_secstruct(work_resid) == 'H' ) {
				helix_length++;
				work_resid++;
			}
			if ( helix_length < 7 ) continue;

			std::string nt = "nterm";

			stmt.bind(1,struct_id);
			stmt.bind(2,i);
			stmt.bind(3,ss);
			stmt.bind(4,nt);
			stmt.bind(5,i-4);
			stmt.bind(6,i-3);
			stmt.bind(7,i-2);
			stmt.bind(8,i-1);
			stmt.bind(9,i+1);
			stmt.bind(10,i+2);
			stmt.bind(11,i+3);
			stmt.bind(12,i+4);
			basic::database::safely_write_to_database(stmt);
		} else {
			continue;
		}
	}


	return 0;
}

} // namespace
} // namespace
