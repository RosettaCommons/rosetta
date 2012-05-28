// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ResidueFeatures.cc
/// @brief  report residue features to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/ResidueFeatures.hh>

//External
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>

// Project Headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

//Basic Headers
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/Constraint.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
//#include <cmath>

namespace protocols{
namespace features{

static basic::Tracer TR("protocols.features.ResidueFeatures");
	
using std::string;
using core::Size;
using core::Real;
using core::pose::Pose;
using core::conformation::Residue;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;


ResidueFeatures::ResidueFeatures() {}

ResidueFeatures::ResidueFeatures( ResidueFeatures const & ) :
	FeaturesReporter()
{}

ResidueFeatures::~ResidueFeatures()
{}

string
ResidueFeatures::type_name() const { return "ResidueFeatures"; }

void
ResidueFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const{
	using namespace basic::database::schema_generator;
	
	Column struct_id("struct_id",DbUUID(), false);
	Column resNum("resNum",DbInteger(), false);
	Column name3("name3",DbText(), false);
	Column res_type("res_type",DbText(), false);
	
	utility::vector1<Column> residues_pkey_cols;
	residues_pkey_cols.push_back(struct_id);
	residues_pkey_cols.push_back(resNum);
	
	Schema residues("residues", PrimaryKey(residues_pkey_cols));
	residues.add_column(struct_id);
	residues.add_column(resNum);
	residues.add_column(name3);
	residues.add_column(res_type);
	residues.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true));
	
	//TODO add constraint resNum > 0
	
	residues.write(db_session);
	
//	if(db_mode == "sqlite3")
//	{
//		return
//			"CREATE TABLE IF NOT EXISTS residues (\n"
//			"	struct_id BLOB,\n"
//			"	resNum INTEGER,\n"
//			"	name3 TEXT,\n"
//			"	res_type TEXT,\n"
//			"	FOREIGN KEY (struct_id)\n"
//			"		REFERENCES structures (struct_id)\n"
//			"		DEFERRABLE INITIALLY DEFERRED,\n"
//			"	CONSTRAINT resNum_is_positive CHECK (resNum >= 1),\n"
//			"	PRIMARY KEY(struct_id, resNum));";
//	}else if(db_mode=="mysql")
//	{
//		return
//			"CREATE TABLE IF NOT EXISTS residues (\n"
//			"	struct_id BINARY(16),\n"
//			"	resNum INTEGER,\n"
//			"	name3 TEXT,\n"
//			"	res_type TEXT,\n"
//			"	FOREIGN KEY (struct_id) REFERENCES structures (struct_id),\n"
//			"	CONSTRAINT resNum_is_positive CHECK (resNum >= 1),\n"
//			"	PRIMARY KEY(struct_id, resNum));";
//	}
}

utility::vector1<std::string>
ResidueFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("StructureFeatures");
	return dependencies;
}

Size
ResidueFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	boost::uuids::uuid const struct_id,
	sessionOP db_session
){
	insert_residue_rows(pose, relevant_residues, struct_id, db_session);
	return 0;
}


void
ResidueFeatures::insert_residue_rows(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	boost::uuids::uuid const struct_id,
	sessionOP db_session
){
	
	std::string statement_string = "INSERT INTO residues (struct_id, resNum, name3, res_type) VALUES (?,?,?,?);";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));

	for(Size resNum=1; resNum <= pose.total_residue(); ++resNum){
		if(!relevant_residues[resNum]) continue;
		Residue res = pose.residue(resNum);

		string const name3( res.name3() );
		string const res_type( res.name() );
		
		//TR << "residues binding - " << to_string(struct_id) << " " << resNum << " " << name3 << " " << res_type << std::endl;

		stmt.bind(1,struct_id);
		stmt.bind(2,resNum);
		stmt.bind(3,name3);
		stmt.bind(4,res_type);
		basic::database::safely_write_to_database(stmt);
	}
}

void
ResidueFeatures::delete_record(
	boost::uuids::uuid struct_id,
	sessionOP db_session) {

	delete_records_from_table("residues", struct_id, db_session);
}



} // namesapce
} // namespace
