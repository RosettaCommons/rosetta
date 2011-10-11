// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/features/PdbDataFeatures.cc
/// @author Sam DeLuca

#include <protocols/features/PdbDataFeatures.hh>

//project headers
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>

// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/database/sql_utils.hh>

//Utility Headers
#include <utility/sql_database/DatabaseSessionManager.hh>

// External Headers
#include <cppdb/frontend.h>



namespace protocols {
namespace features {

PdbDataFeatures::PdbDataFeatures()
{

}

PdbDataFeatures::PdbDataFeatures(PdbDataFeatures const & src)
{

}

PdbDataFeatures::~PdbDataFeatures()
{

}

std::string PdbDataFeatures::type_name() const
{
	return "PdbDataFeatures";
}

std::string PdbDataFeatures::schema() const
{
	std::string db_mode(basic::options::option[basic::options::OptionKeys::inout::database_mode]);
	if(db_mode == "sqlite3")
	{
		return
			"CREATE TABLE IF NOT EXISTS residue_pdb_identification (\n"
			"	struct_id INTEGER,\n"
			"	residue_number INTEGER,\n"
			"	chain_id TEXT,\n"
			"	insertion_code TEXT,\n"
			"	pdb_residue_number INTEGER,\n"
			"	FOREIGN KEY (struct_id)\n"
			"	REFERENCES structures(struct_id)"
			"	DEFERRABLE INITIALLY DEFERRED);";
	}else if(db_mode=="mysql")
	{
		return
			"CREATE TABLE IF NOT EXISTS residue_pdb_identification (\n"
			"	struct_id INTEGER,\n"
			"	residue_number INTEGER,\n"
			"	chain_id TEXT,\n"
			"	insertion_code TEXT,\n"
			"	pdb_residue_number INTEGER,\n"
			"	FOREIGN KEY (struct_id)\n"
			"	REFERENCES structures(struct_id)"
			"	DEFERRABLE INITIALLY DEFERRED);";
	}else
	{
		return "";
	}
}

core::Size PdbDataFeatures::report_features(
	core::pose::Pose const & pose,
	utility::vector1<bool> const & relevant_residues,
	core::Size struct_id,
	utility::sql_database::sessionOP db_session )
{
	insert_pdb_info_rows(struct_id,db_session,pose);
	return 0;
}

void PdbDataFeatures::delete_record(
	core::Size struct_id,
	utility::sql_database::sessionOP db_session)
{
	cppdb::statement stmt = (*db_session) << "DELETE FROM residue_pdb_identification WHERE struct_id == ?;\n" <<struct_id;
	stmt.exec();
}

void PdbDataFeatures::load_into_pose(
	utility::sql_database::sessionOP db_session,
	core::Size struct_id,
	core::pose::Pose & pose)
{
	load_pdb_info(db_session,struct_id,pose);
}

void PdbDataFeatures::load_pdb_info(
	utility::sql_database::sessionOP db_session,
	core::Size struct_id,
	core::pose::Pose & pose)
{

	utility::vector1<core::Size> pdb_numbers;
	utility::vector1<char> pdb_chains;
	utility::vector1<char> insertion_codes;
	cppdb::statement stmt = (*db_session) <<
		"SELECT\n"
		"	residue_number,\n"
		"	chain_id,\n"
		"	insertion_code,\n"
		"	pdb_residue_number\n"
		"FROM\n"
		"	residue_pdb_identification\n"
		"WHERE\n"
		"	residue_pdb_identification.struct_id=?;" <<struct_id;

	cppdb::result res(basic::database::safely_read_from_database(stmt));

	while(res.next()) {
		core::Size residue_number;
		//cppdb doesn't do char's
		std::string chain_id;
		std::string insertion_code;
		core::Size pdb_residue_number;

		res >> residue_number >> chain_id >> insertion_code >> pdb_residue_number;

		pdb_numbers.push_back(residue_number);
		pdb_chains.push_back(chain_id[0]);
		insertion_codes.push_back(insertion_code[0]);
	}

	core::pose::PDBInfoOP pdb_info( new core::pose::PDBInfo( pose.total_residue() ) );

	pdb_info->set_numbering(pdb_numbers);
	pdb_info->set_chains(pdb_chains);
	pdb_info->set_icodes(insertion_codes);

	pose.pdb_info(pdb_info);

}

void PdbDataFeatures::insert_pdb_info_rows(core::Size struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose)
{
	core::Size res_num(pose.n_residue());
	for(core::Size index =1 ; index <= res_num; ++index)
	{
		std::string chain_id(& pose.pdb_info()->chain(index),1);
		std::string insertion_code(&pose.pdb_info()->icode(index),1);
		core::Size pdb_residue_number = pose.pdb_info()->number(index);

		cppdb::statement stmt = (*db_session)
			<< "INSERT INTO residue_pdb_identification VALUES (?,?,?,?,?);"
			<< struct_id
			<< index
			<< chain_id
			<< insertion_code
			<< pdb_residue_number;
		basic::database::safely_read_from_database(stmt);
	}
}

}
}


