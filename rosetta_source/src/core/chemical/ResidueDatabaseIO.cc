// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/chemical/residue_dbio.cc
/// @author Sam DeLuca
/// @author Matt O'Meara

// Project Headers
#include <basic/Tracer.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueDatabaseIO.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/kinematics/Stub.hh>
#include <core/types.hh>
#include <basic/database/sql_utils.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>

#include <numeric/xyz.functions.hh>

#include <utility/string_util.hh>
// Boost Headers
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

namespace core {
namespace chemical {

ResidueDatabaseIO::ResidueDatabaseIO() :version_(0.2)
{

}

ResidueDatabaseIO::~ResidueDatabaseIO()
{

}


std::string ResidueDatabaseIO::schema() const
{
	// NOTE: To support building feature databases in parallel, the
	// ResidueTypeSet and ResidueType objects must be identified by
	// their names rather then assigning them a unique id.

	std::string db_mode(basic::options::option[basic::options::OptionKeys::inout::database_mode]);
	if(db_mode == "sqlite3")
	{
		return
			"CREATE TABLE IF NOT EXISTS residue_type (\n"
			"	residue_type_set_name TEXT,\n"
			"	version TEXT,\n"
			"	name TEXT,\n"
			"	name3 TEXT,\n"
			"	name1 TEXT,\n"
			"	aa TEXT,\n"
			"	lower_connect INTEGER,\n"
			"	upper_connect INTEGER,\n"
			"	nbr_atom INTEGER,\n"
			"	nbr_radius REAL,\n"
			"	rotamer_library TEXT,\n"
			"	FOREIGN KEY(residue_type_set_name)\n"
			"		REFERENCES residue_type_set(name)\n"
			"		DEFERRABLE INITIALLY DEFERRED,\n"
			"	PRIMARY KEY(residue_type_set_name, name));\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS residue_type_atom (\n"
			"	residue_type_set_name TEXT,\n"
			"	residue_type_name TEXT,\n"
			"	atom_index INTEGER,\n"
			"	atom_name TEXT,\n"
			"	atom_type_name TEXT,\n"
			"	mm_atom_type_name TEXT,\n"
			"	charge REAL,\n"
			"	is_backbone INTEGER,\n"
			"	FOREIGN KEY(residue_type_set_name, residue_type_name)\n"
			"		REFERENCES residue_type(residue_type_set_name, name)\n"
			"		DEFERRABLE INITIALLY DEFERRED,\n"
			"	PRIMARY KEY(residue_type_set_name, residue_type_name, atom_index));\n"
			"\n"
			"CREATE UNIQUE INDEX IF NOT EXISTS\n"
			"	residue_type_atom_residue_type_set_name_residue_type_name_atom_name ON\n"
			"	residue_type_atom ( residue_type_set_name, residue_type_name, atom_name );\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS residue_type_bond (\n"
			"	residue_type_set_name TEXT,\n"
			"	residue_type_name TEXT,\n"
			"	atom1 INTEGER,\n"
			"	atom2 INTEGER,\n"
			"	bond_type INTEGER,\n"
			"	FOREIGN KEY(residue_type_set_name, residue_type_name)\n"
			"		REFERENCES residue_type(residue_type_set_name, name)\n"
			"		DEFERRABLE INITIALLY DEFERRED,\n"
			"	PRIMARY KEY(residue_type_set_name, residue_type_name, atom1, atom2));\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS residue_type_cut_bond (\n"
			"	residue_type_set_name TEXT,\n"
			"	residue_type_name TEXT,\n"
			"	atom1 INTEGER,\n"
			"	atom2 INTEGER,\n"
			"	FOREIGN KEY(residue_type_set_name, residue_type_name)\n"
			"		REFERENCES residue_type(residue_type_set_name, name)\n"
			"		DEFERRABLE INITIALLY DEFERRED,\n"
			"	PRIMARY KEY(residue_type_set_name, residue_type_name, atom1, atom2));\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS residue_type_chi (\n"
			"	residue_type_set_name TEXT,\n"
			"	residue_type_name TEXT,\n"
			"	chino INTEGER,\n"
			"	atom1 TEXT,\n"
			"	atom2 TEXT,\n"
			"	atom3 TEXT,\n"
			"	atom4 TEXT,\n"
			"	FOREIGN KEY(residue_type_set_name, residue_type_name)\n"
			"		REFERENCES residue_type(residue_type_set_name, name)\n"
			"		DEFERRABLE INITIALLY DEFERRED,\n"
			"	PRIMARY KEY(residue_type_set_name, residue_type_name, atom1, atom2, atom3, atom4));\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS residue_type_chi_rotamer (\n"
			"	residue_type_set_name TEXT,\n"
			"	residue_type_name TEXT,\n"
			"	chino INTEGER,\n"
			"	mean REAL,\n"
			"	sdev REAL,\n"
			"	FOREIGN KEY(residue_type_set_name, residue_type_name)\n"
			"		REFERENCES residue_type(residue_type_set_name, name)\n"
			"		DEFERRABLE INITIALLY DEFERRED,\n"
			"	PRIMARY KEY(residue_type_set_name, residue_type_name, chino));\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS residue_type_proton_chi (\n"
			"	residue_type_set_name TEXT,\n"
			"	residue_type_name TEXT,\n"
			"	chino INTEGER,\n"
			"	sample REAL,\n"
			"	is_extra BOOL,\n"
			"	FOREIGN KEY(residue_type_set_name, residue_type_name)\n"
			"		REFERENCES residue_type(residue_type_set_name, name)\n"
			"		DEFERRABLE INITIALLY DEFERRED,\n"
			"	PRIMARY KEY(residue_type_set_name, residue_type_name, chino, sample));\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS residue_type_property (\n"
			"	residue_type_set_name TEXT,\n"
			"	residue_type_name TEXT,\n"
			"	property TEXT,\n"
			"	FOREIGN KEY(residue_type_set_name, residue_type_name)\n"
			"		REFERENCES residue_type(residue_type_set_name, name)\n"
			"		DEFERRABLE INITIALLY DEFERRED,\n"
			"	PRIMARY KEY(residue_type_set_name, residue_type_name, property));\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS residue_type_variant_type (\n"
			"	residue_type_set_name TEXT,\n"
			"	residue_type_name TEXT,\n"
			"	variant_type TEXT,\n"
			"	FOREIGN KEY(residue_type_set_name, residue_type_name)\n"
			"		REFERENCES residue_type(residue_type_set_name, name)\n"
			"		DEFERRABLE INITIALLY DEFERRED,\n"
			"	PRIMARY KEY(residue_type_set_name, residue_type_name, variant_type));\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS residue_type_icoor (\n"
			"	residue_type_set_name TEXT,\n"
			"	residue_type_name TEXT,\n"
			"	icoor_sequence INTEGER,\n"
			"	child_atom TEXT,\n"
			"	phi REAL,\n"
			"	theta REAL,\n"
			"	distance REAL,\n"
			"	parent_atom TEXT,\n"
			"	angle_atom TEXT,\n"
			"	torsion_atom TEXT,\n"
			"	FOREIGN KEY(residue_type_set_name,residue_type_name)\n"
			"		REFERENCES residue_type(residue_type_set_name,name)\n"
			"		DEFERRABLE INITIALLY DEFERRED,\n"
			"	PRIMARY KEY(residue_type_set_name,residue_type_name,child_atom));\n";

	}else if (db_mode == "mysql")
	{
		return
			"CREATE TABLE IF NOT EXISTS residue_type (\n"
			"	residue_type_set_name VARCHAR(255),\n"
			"	version TEXT,\n"
			"	name VARCHAR(255),\n"
			"	name3 TEXT,\n"
			"	name1 TEXT,\n"
			"	aa TEXT,\n"
			"	lower_connect INTEGER,\n"
			"	upper_connect INTEGER,\n"
			"	nbr_atom INTEGER,\n"
			"	nbr_radius REAL,\n"
			"	rotamer_library TEXT,\n"
			"	PRIMARY KEY(residue_type_set_name, name));\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS residue_type_atom (\n"
			"	residue_type_set_name VARCHAR(255),\n"
			"	residue_type_name VARCHAR(255),\n"
			"	atom_index INTEGER,\n"
			"	atom_name TEXT,\n"
			"	atom_type_name TEXT,\n"
			"	mm_atom_type_name TEXT,\n"
			"	charge REAL,\n"
			"	is_backbone INTEGER,\n"
			"	FOREIGN KEY(residue_type_set_name, residue_type_name) REFERENCES residue_type(residue_type_set_name, name),\n"
			"	PRIMARY KEY(residue_type_set_name, residue_type_name, atom_index));\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS residue_type_bond (\n"
			"	residue_type_set_name VARCHAR(255),\n"
			"	residue_type_name VARCHAR(255),\n"
			"	atom1 INTEGER,\n"
			"	atom2 INTEGER,\n"
			"	bond_type INTEGER,\n"
			"	FOREIGN KEY(residue_type_set_name, residue_type_name) REFERENCES residue_type(residue_type_set_name, name),\n"
			"	PRIMARY KEY(residue_type_set_name, residue_type_name, atom1, atom2));\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS residue_type_cut_bond (\n"
			"	residue_type_set_name VARCHAR(255),\n"
			"	residue_type_name VARCHAR(255),\n"
			"	atom1 INTEGER,\n"
			"	atom2 INTEGER,\n"
			"	FOREIGN KEY(residue_type_set_name, residue_type_name) REFERENCES residue_type(residue_type_set_name, name),\n"
			"	PRIMARY KEY(residue_type_set_name, residue_type_name, atom1, atom2));\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS residue_type_chi (\n"
			"	residue_type_set_name VARCHAR(255),\n"
			"	residue_type_name VARCHAR(255),\n"
			"	chino INTEGER,\n"
			"	atom1 VARCHAR(8),\n"
			"	atom2 VARCHAR(8),\n"
			"	atom3 VARCHAR(8),\n"
			"	atom4 VARCHAR(8),\n"
			"	FOREIGN KEY(residue_type_set_name, residue_type_name) REFERENCES residue_type(residue_type_set_name, name),\n"
			"	PRIMARY KEY(residue_type_set_name, residue_type_name, atom1, atom2, atom3, atom4));\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS residue_type_chi_rotamer (\n"
			"	residue_type_set_name VARCHAR(255),\n"
			"	residue_type_name VARCHAR(255),\n"
			"	chino INTEGER,\n"
			"	mean REAL,\n"
			"	sdev REAL,\n"
			"	FOREIGN KEY(residue_type_set_name, residue_type_name) REFERENCES residue_type(residue_type_set_name, name),\n"
			"	PRIMARY KEY(residue_type_set_name, residue_type_name, chino));\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS residue_type_proton_chi (\n"
			"	residue_type_set_name VARCHAR(255),\n"
			"	residue_type_name VARCHAR(255),\n"
			"	chino INTEGER,\n"
			"	sample REAL,\n"
			"	is_extra BOOL,\n"
			"	FOREIGN KEY(residue_type_set_name, residue_type_name) REFERENCES residue_type(residue_type_set_name, name),\n"
			"	PRIMARY KEY(residue_type_set_name, residue_type_name, chino, sample));\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS residue_type_property (\n"
			"	residue_type_set_name VARCHAR(255),\n"
			"	residue_type_name VARCHAR(255),\n"
			"	property VARCHAR(255),\n"
			"	FOREIGN KEY(residue_type_set_name, residue_type_name) REFERENCES residue_type(residue_type_set_name, name),\n"
			"	PRIMARY KEY(residue_type_set_name, residue_type_name, property));\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS residue_type_variant_type (\n"
			"	residue_type_set_name VARCHAR(255),\n"
			"	residue_type_name VARCHAR(255),\n"
			"	variant_type VARCHAR(255),\n"
			"	FOREIGN KEY(residue_type_set_name, residue_type_name) REFERENCES residue_type(residue_type_set_name, name),\n"
			"	PRIMARY KEY(residue_type_set_name, residue_type_name, variant_type));\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS residue_type_icoor (\n"
			"	residue_type_set_name VARCHAR(255),\n"
			"	residue_type_name VARCHAR(255),\n"
			"	icoor_sequence INTEGER,\n"
			"	child_atom VARCHAR(8),\n"
			"	phi REAL,\n"
			"	theta REAL,\n"
			"	distance REAL,\n"
			"	parent_atom VARCHAR(8),\n"
			"	angle_atom VARCHAR(8),\n"
			"	torsion_atom VARCHAR(8),\n"
			"	FOREIGN KEY(residue_type_set_name,residue_type_name) REFERENCES residue_type(residue_type_set_name,name),\n"
			"	PRIMARY KEY(residue_type_set_name,residue_type_name,child_atom));\n";

	}else
	{
		return "";
	}

}

void ResidueDatabaseIO::initialize(utility::sql_database::sessionOP db_session)
{
	basic::database::write_schema_to_database(schema(),db_session);
}

core::Real ResidueDatabaseIO::get_version()
{
	return version_;
}


void ResidueDatabaseIO::write_residuetype_to_database(
	std::string const & residue_type_set_name,
	core::chemical::ResidueType const & res_type,
	utility::sql_database::sessionOP db_session)
{
	std::string status_string = "SELECT * FROM residue_type WHERE residue_type_set_name = ? AND name = ?;";
	cppdb::statement status_statement(basic::database::safely_prepare_statement(status_string,db_session));
	status_statement.bind(1,residue_type_set_name);
	status_statement.bind(2,res_type.name());
	cppdb::result res(basic::database::safely_read_from_database(status_statement));

	if(res.next()) return;

	report_residue_type(residue_type_set_name, res_type, db_session);
	report_residue_type_atom(residue_type_set_name, res_type, db_session);
	report_residue_type_bond(residue_type_set_name, res_type, db_session);
	report_residue_type_cut_bond(residue_type_set_name, res_type, db_session);
	report_residue_type_chi(residue_type_set_name, res_type, db_session);
	report_residue_type_chi_rotamer(residue_type_set_name, res_type, db_session);
	report_residue_type_proton_chi(residue_type_set_name, res_type, db_session);
	report_residue_type_properties(residue_type_set_name, res_type, db_session);
	report_residue_type_variant(residue_type_set_name, res_type, db_session);
	report_residue_type_icoor(residue_type_set_name, res_type, db_session);
}

core::chemical::ResidueTypeOP ResidueDatabaseIO::read_residuetype_from_database(
	chemical::AtomTypeSetCAP atom_types,
	chemical::ElementSetCAP elements,
	chemical::MMAtomTypeSetCAP mm_atom_types,
	chemical::orbitals::OrbitalTypeSetCAP orbital_atom_types,
	std::string const & residue_type_set_name,
	std::string const & residue_type_name,
	utility::sql_database::sessionOP db_session
)
{

	core::chemical::ResidueTypeOP res_type(new core::chemical::ResidueType(atom_types,elements,mm_atom_types,orbital_atom_types));


	read_residue_type_atom(residue_type_set_name,residue_type_name,*res_type,db_session);
	read_residue_type(residue_type_set_name,residue_type_name,*res_type,db_session);
	read_residue_type_bond(residue_type_set_name,residue_type_name,*res_type,db_session);
	read_residue_type_cut_bond(residue_type_set_name,residue_type_name,*res_type,db_session);
	read_residue_type_chi(residue_type_set_name,residue_type_name,*res_type,db_session);
	read_residue_type_chi_rotamer(residue_type_set_name,residue_type_name,*res_type,db_session);
	read_residue_type_properties(residue_type_set_name,residue_type_name,*res_type,db_session);
	read_residue_type_variant(residue_type_set_name,residue_type_name,*res_type,db_session);
	read_residue_type_icoor(residue_type_set_name,residue_type_name,*res_type,db_session);
	res_type->finalize();

	return res_type;
}

std::string ResidueDatabaseIO::get_atom_name_from_database_atom_index(
	std::string residue_name,
	core::Size atom_index,
	utility::sql_database::sessionOP db_session
)  {

	std::pair<std::string,core::Size> atom_id = std::make_pair(residue_name,atom_index);

	if(atom_name_id_cache_.find(atom_id) != atom_name_id_cache_.end())
	{
		return atom_name_id_cache_[atom_id];
	}

	std::string atom_query =
		"SELECT\n"
		"	atom_index,\n"
		"	atom_name\n"
		"FROM residue_type_atom WHERE residue_type_name=?;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(atom_query,db_session));
	stmt.bind(1,residue_name);
	cppdb::result res(basic::database::safely_read_from_database(stmt));

	while(res.next())
	{
		std::string atom_name;
		core::Size atom_index;

		res >> atom_index >> atom_name;

		std::pair<std::string,core::Size> new_atom_id = std::make_pair(residue_name,atom_index);
		atom_name_id_cache_.insert(std::make_pair(new_atom_id,atom_name));

	}

	return atom_name_id_cache_[atom_id];


}

void ResidueDatabaseIO::report_residue_type(
	std::string const & residue_type_set_name,
	ResidueType const & res_type,
	utility::sql_database::sessionOP db_session
)  {

	std::stringstream name1; name1 << res_type.name1();

	int lower_terminus(-1), upper_terminus(-1);
	if(res_type.is_polymer()){
		if(!res_type.is_lower_terminus()) lower_terminus = res_type.lower_connect_atom();
		if(!res_type.is_upper_terminus()) upper_terminus = res_type.upper_connect_atom();
	}

	std::string statement_string = "INSERT INTO residue_type VALUES (?,?,?,?,?,?,?,?,?,?,?);";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));

	stmt.bind(1,residue_type_set_name);
	stmt.bind(2,version_);
	stmt.bind(3,res_type.name());
	stmt.bind(4,res_type.name3());
	stmt.bind(5,name1.str());
	stmt.bind(6,res_type.aa());
	stmt.bind(7,lower_terminus);
	stmt.bind(8,upper_terminus);
	stmt.bind(9,res_type.nbr_atom());
	stmt.bind(10,res_type.nbr_radius());
	stmt.bind(11,res_type.get_RotamerLibraryName());
	basic::database::safely_write_to_database(stmt);

}

/// @detail this needs to get called after read_residue_type_atom
void ResidueDatabaseIO::read_residue_type(
		std::string const & residue_type_set_name,
		std::string const & residue_type_name,
		core::chemical::ResidueType & res_type,
		utility::sql_database::sessionOP db_session)
{
	std::string residue_type_statement =
		"SELECT\n"
		"	version,\n"
		"	name3,\n"
		"	name1,\n"
		"	aa,\n"
		"	lower_connect,\n"
		"	upper_connect,\n"
		"	nbr_atom,\n"
		"	nbr_radius,\n"
		"	rotamer_library\n"
		"FROM residue_type\n"
		"WHERE residue_type_set_name = ? AND name = ?;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(residue_type_statement,db_session));
	stmt.bind(1,residue_type_set_name);
	stmt.bind(2,residue_type_name);
	cppdb::result res(basic::database::safely_read_from_database(stmt));

	if(!res.next())
	{
		utility_exit_with_message("could not find residue "+residue_type_name+" in "+residue_type_set_name);
	}

	core::Real version;
	std::string name3, name1;
	core::Size aa;
	int lower_connect, upper_connect;
	core::Size nbr_atom;
	std::string rotamer_library;
	core::Real nbr_radius;

	res >>
		version >>
		name3 >>
		name1 >>
		aa >>
		lower_connect >>
		upper_connect >>
		nbr_atom >>
		nbr_radius >>
		rotamer_library;
	if(version != version_)
	{
		utility_exit_with_message("Version mismatch between Residue Database and Executable");
	}

	res_type.name(residue_type_name);
	res_type.name3(name3);
	res_type.name1(name1[0]);
	res_type.aa(name_from_aa(static_cast<AA>(aa)));

	if(lower_connect > 0 )
	{
		res_type.set_lower_connect_atom(get_atom_name_from_database_atom_index(residue_type_name,lower_connect,db_session));
	}

	if(upper_connect > 0 )
	{
		res_type.set_upper_connect_atom(get_atom_name_from_database_atom_index(residue_type_name,lower_connect,db_session));
	}

	res_type.nbr_atom(get_atom_name_from_database_atom_index(residue_type_name,nbr_atom,db_session));
	res_type.nbr_radius(nbr_radius);

	if(rotamer_library != "")
	{
		utility::file::FileName rot_file( rotamer_library );
		res_type.set_RotamerLibraryName(rotamer_library);
	}

}

utility::vector1<std::string> ResidueDatabaseIO::get_all_residues_in_database(utility::sql_database::sessionOP db_session)
{

	utility::vector1<std::string> residue_names;

	std::string residue_name_statement =
		"SELECT\n"
		"	name\n"
		"FROM residue_type;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(residue_name_statement,db_session));
	cppdb::result res(basic::database::safely_read_from_database(stmt));
	while(res.next())
	{
		std::string name;
		res >> name;
		residue_names.push_back(name);
	}
	return residue_names;
}

void
ResidueDatabaseIO::report_residue_type_atom(
	std::string const & residue_type_set_name,
	ResidueType const & res_type,
	utility::sql_database::sessionOP db_session
)  {

	std::string statement_string = "INSERT INTO residue_type_atom VALUES (?,?,?,?,?,?,?,?);";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	// AtomTypeSet?
	for(Size i=1; i <= res_type.natoms(); ++i){

		stmt.bind(1,residue_type_set_name);
		stmt.bind(2,res_type.name());
		stmt.bind(3,i);
		stmt.bind(4,res_type.atom_name(i));
		stmt.bind(5,res_type.atom_type(i).atom_type_name());
		stmt.bind(6,res_type.mm_atom_name(i));
		stmt.bind(7,res_type.atomic_charge(i));
		stmt.bind(8,res_type.atom_is_backbone(i));
		basic::database::safely_write_to_database(stmt);

	}
}

void
ResidueDatabaseIO::read_residue_type_atom(
	std::string const & residue_type_set_name,
	std::string const & residue_type_name,
	core::chemical::ResidueType & res_type,
	utility::sql_database::sessionOP db_session)
{
	std::string residue_atoms_statement =
		"SELECT\n"
		"atom_index,\n"
		"atom_name,\n"
		"atom_type_name,\n"
		"mm_atom_type_name,\n"
		"charge\n"
		"FROM residue_type_atom\n"
		"WHERE residue_type_set_name = ? AND residue_type_name = ?\n"
		"ORDER BY atom_index;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(residue_atoms_statement,db_session));
	stmt.bind(1,residue_type_set_name);
	stmt.bind(2,residue_type_name);
	cppdb::result res(basic::database::safely_read_from_database(stmt));

	while(res.next())
	{
		core::Size atom_index;
		std::string atom_name, atom_type_name,mm_atom_type_name;
		core::Real charge;
		core::Size is_backbone;

		res >> atom_index >> atom_name >> atom_type_name >> mm_atom_type_name >> charge;
		res_type.add_atom(atom_name,atom_type_name,mm_atom_type_name,charge);
	}

}

void
ResidueDatabaseIO::report_residue_type_bond(
	std::string const & residue_type_set_name,
	ResidueType const & res_type,
	utility::sql_database::sessionOP db_session
)  {

	std::string statement_string = "INSERT INTO residue_type_bond VALUES (?,?,?,?,?);";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	for(Size atm=1; atm <= res_type.natoms(); ++atm){
		AtomIndices const & neighbors(res_type.bonded_neighbor(atm));
		utility::vector1<BondName> const & types(res_type.bonded_neighbor_types(atm));
		for(Size nbr=1; nbr <= neighbors.size(); ++nbr){
			if(atm >= neighbors[nbr]) continue;

			stmt.bind(1,residue_type_set_name);
			stmt.bind(2,res_type.name());
			stmt.bind(3,atm);
			stmt.bind(4,neighbors[nbr]);
			stmt.bind(5,types[nbr]);
			basic::database::safely_write_to_database(stmt);
		}
	}
}

void
ResidueDatabaseIO::read_residue_type_bond(
	std::string const & residue_type_set_name,
	std::string const & residue_type_name,
	core::chemical::ResidueType & res_type,
	utility::sql_database::sessionOP db_session
)
{
	std::string bond_statement_string =
		"SELECT\n"
		"	atom1,\n"
		"	atom2,\n"
		"	bond_type\n"
		"FROM residue_type_bond\n"
		"WHERE residue_type_set_name = ? AND residue_type_name = ?;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(bond_statement_string,db_session));
	stmt.bind(1,residue_type_set_name);
	stmt.bind(2,residue_type_name);
	cppdb::result res(basic::database::safely_read_from_database(stmt));
	while(res.next())
	{
		core::Size atom1, atom2, bond_type;
		res >> atom1 >> atom2 >> bond_type;
		std::string atom1_name = get_atom_name_from_database_atom_index(residue_type_name,atom1,db_session);
		std::string atom2_name = get_atom_name_from_database_atom_index(residue_type_name,atom2,db_session);
		res_type.add_bond(atom1_name,atom2_name,static_cast<core::chemical::BondName>(bond_type));
	}


}

void
ResidueDatabaseIO::report_residue_type_cut_bond(
	std::string const & residue_type_set_name,
	ResidueType const & res_type,
	utility::sql_database::sessionOP db_session
)  {

	std::string statement_string = "INSERT INTO residue_type_cut_bond VALUES (?,?,?,?);";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));

	for(Size i=1; i <= res_type.natoms(); ++i){
		foreach(core::Size const j, res_type.cut_bond_neighbor(i)){
			if(i>=j) continue;


			stmt.bind(1,residue_type_set_name);
			stmt.bind(2,res_type.name());
			stmt.bind(3,i);
			stmt.bind(4,j);
			basic::database::safely_write_to_database(stmt);
		}
	}
}

void
ResidueDatabaseIO::read_residue_type_cut_bond(
	std::string const & residue_type_set_name,
	std::string const & residue_type_name,
	core::chemical::ResidueType & res_type,
	utility::sql_database::sessionOP db_session
)
{
	std::string cut_bond_statement_string =
		"SELECT\n"
		"	atom1,\n"
		"	atom2\n"
		"FROM residue_type_cut_bond\n"
		"WHERE residue_type_set_name = ? AND residue_type_name = ?;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(cut_bond_statement_string,db_session));
	stmt.bind(1,residue_type_set_name);
	stmt.bind(2,residue_type_name);
	cppdb::result res(basic::database::safely_read_from_database(stmt));
	while(res.next())
	{
		core::Size atom1, atom2;
		res >> atom1 >> atom2;
		std::string atom1_name = get_atom_name_from_database_atom_index(residue_type_name,atom1,db_session);
		std::string atom2_name = get_atom_name_from_database_atom_index(residue_type_name,atom2,db_session);
		res_type.add_cut_bond(atom1_name,atom2_name);
	}
}


void
ResidueDatabaseIO::report_residue_type_chi(
	std::string const & residue_type_set_name,
	ResidueType const & res_type,
	utility::sql_database::sessionOP db_session
)  {

	std::string statement_string = "INSERT INTO residue_type_chi VALUES (?,?,?,?,?,?,?);";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));

	for(Size i=1; i <= res_type.nchi(); ++i){
		AtomIndices const & chi_atoms(res_type.chi_atoms(i));

		stmt.bind(1,residue_type_set_name);
		stmt.bind(2,res_type.name());
		stmt.bind(3,i);
		stmt.bind(4,chi_atoms[1]);
		stmt.bind(5,chi_atoms[2]);
		stmt.bind(6,chi_atoms[3]);
		stmt.bind(7,chi_atoms[4]);
		basic::database::safely_write_to_database(stmt);
	}

}


void
ResidueDatabaseIO::read_residue_type_chi(
	std::string const & residue_type_set_name,
	std::string const & residue_type_name,
	core::chemical::ResidueType & res_type,
	utility::sql_database::sessionOP db_session
)
{
	std::string chi_statement_string =
		"SELECT\n"
		"	chino,\n"
		"	atom1,\n"
		"	atom2,\n"
		"	atom3,\n"
		"	atom4\n"
		"FROM residue_type_chi\n"
		"WHERE residue_type_set_name = ? AND residue_type_name = ?;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(chi_statement_string,db_session));
	stmt.bind(1,residue_type_set_name);
	stmt.bind(2,residue_type_name);
	cppdb::result res(basic::database::safely_read_from_database(stmt));
	while(res.next())
	{
		core::Size chino;
		core::Size atom1, atom2, atom3, atom4;
		res >> chino >> atom1 >> atom2 >> atom3 >> atom4;
		std::string atom1_name = get_atom_name_from_database_atom_index(residue_type_name,atom1,db_session);
		std::string atom2_name = get_atom_name_from_database_atom_index(residue_type_name,atom2,db_session);
		std::string atom3_name = get_atom_name_from_database_atom_index(residue_type_name,atom3,db_session);
		std::string atom4_name = get_atom_name_from_database_atom_index(residue_type_name,atom4,db_session);
		res_type.add_chi(chino,atom1_name,atom2_name,atom3_name,atom4_name);
	}
}


void
ResidueDatabaseIO::report_residue_type_chi_rotamer(
	std::string const & residue_type_set_name,
	ResidueType const & res_type,
	utility::sql_database::sessionOP db_session
)  {

	std::string statement_string = "INSERT INTO residue_type_chi_rotamer VALUES (?,?,?,?,?);";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	for(Size chi=1; chi <= res_type.nchi(); ++chi){
		std::pair<Real, Real> mean_sdev;
		foreach(mean_sdev, res_type.chi_rotamers(chi)){

			stmt.bind(1,residue_type_set_name);
			stmt.bind(2,res_type.name());
			stmt.bind(3,chi);
			stmt.bind(4,mean_sdev.first);
			stmt.bind(5,mean_sdev.second);
			basic::database::safely_write_to_database(stmt);
		}
	}
}

void
ResidueDatabaseIO::read_residue_type_chi_rotamer(
	std::string const & residue_type_set_name,
	std::string const & residue_type_name,
	core::chemical::ResidueType & res_type,
	utility::sql_database::sessionOP db_session)
{
	std::string chi_rotamer_statement_string =
		"SELECT\n"
		"	chino,\n"
		"	mean,\n"
		"	sdev\n"
		"FROM residue_type_chi_rotamer\n"
		"WHERE residue_type_set_name = ? AND residue_type_name = ?;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(chi_rotamer_statement_string,db_session));
	stmt.bind(1,residue_type_set_name);
	stmt.bind(2,residue_type_name);
	cppdb::result res(basic::database::safely_read_from_database(stmt));
	while(res.next())
	{
		core::Size chino;
		core::Real mean, stdev;

		res >> chino >> mean >> stdev;
		res_type.add_chi_rotamer(chino,mean,stdev);
	}
}

void
ResidueDatabaseIO::report_residue_type_proton_chi(
	std::string const & residue_type_set_name,
	ResidueType const & res_type,
	utility::sql_database::sessionOP db_session
)  {

	std::string statement_string = "INSERT INTO residue_type_proton_chi VALUES (?,?,?,?,?);";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	for(Size proton_chi=1; proton_chi <= res_type.n_proton_chi(); ++proton_chi){
		Size const chi(res_type.proton_chi_2_chi(proton_chi));
		foreach(Real const sample, res_type.proton_chi_samples(proton_chi)){

			stmt.bind(1,residue_type_set_name);
			stmt.bind(2,res_type.name());
			stmt.bind(3,chi);
			stmt.bind(4,sample);
			stmt.bind(5,false);
			basic::database::safely_write_to_database(stmt);

			foreach(Real const extra_sample, res_type.proton_chi_extra_samples(proton_chi)){

				stmt.bind(1,residue_type_set_name);
				stmt.bind(2,res_type.name());
				stmt.bind(3,chi);
				stmt.bind(4,sample - extra_sample);
				stmt.bind(5,true);
				basic::database::safely_write_to_database(stmt);

				stmt.bind(1,residue_type_set_name);
				stmt.bind(2,res_type.name());
				stmt.bind(3,chi);
				stmt.bind(4,sample + extra_sample);
				stmt.bind(5,true);
				basic::database::safely_write_to_database(stmt);

			}
		}
	}
}

void
ResidueDatabaseIO::read_residue_type_proton_chi(
	std::string const & residue_type_set_name,
	std::string const & residue_type_name,
	core::chemical::ResidueType & res_type,
	utility::sql_database::sessionOP db_session)
{
	std::string proton_chi_statement_string =
		"SELECT\n"
		"	chino,\n"
		"	sample,\n"
		"	is_extra\n"
		"FROM residue_type_proton_chi\n"
		"WHERE residue_type_set_name = ? AND residue_type_name = ?;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(proton_chi_statement_string,db_session));
	stmt.bind(1,residue_type_set_name);
	stmt.bind(2,residue_type_name);
	cppdb::result res(basic::database::safely_read_from_database(stmt));

	std::map<core::Size,utility::vector1< core::Real > > samples;
	std::map<core::Size,utility::vector1< core::Real > > extra_samples;

	while(res.next())
	{
		core::Size chino, is_extra;
		core::Real sample;

		res >> chino >> sample >> is_extra;
		if(is_extra)
		{
			if(extra_samples.find(chino) == extra_samples.end())
			{
				utility::vector1< core::Real> extra_sample_vect;
				extra_sample_vect.push_back(sample);
				extra_samples.insert(std::make_pair(chino,extra_sample_vect));
			}else
			{
				extra_samples[chino].push_back(sample);
			}
		}else
		{
			//we've never seen this chino, add it to the map
			if(samples.find(chino) == samples.end())
			{
				utility::vector1< core::Real> sample_vect;
				sample_vect.push_back(sample);
				samples.insert(std::make_pair(chino,sample_vect));
			}else
			{
				samples[chino].push_back(sample);
			}
		}
	}
}

void
ResidueDatabaseIO::report_residue_type_properties(
	std::string const & residue_type_set_name,
	ResidueType const & res_type,
	utility::sql_database::sessionOP db_session
)  {
	foreach(std::string const & property, res_type.properties()){

		cppdb::statement stmt = (*db_session)
					<< "INSERT INTO residue_type_property VALUES (?,?,?);"
					<< residue_type_set_name
					<< res_type.name()
					<< property;
		basic::database::safely_write_to_database(stmt);

	}
}

void
ResidueDatabaseIO::read_residue_type_properties(
	std::string const & residue_type_set_name,
	std::string const & residue_type_name,
	core::chemical::ResidueType & res_type,
	utility::sql_database::sessionOP db_session)
{
	std::string residue_property_statement_string =
		"SELECT\n"
		"	property\n"
		"FROM residue_type_property\n"
		"WHERE residue_type_set_name = ? AND residue_type_name = ?;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(residue_property_statement_string,db_session));
	stmt.bind(1,residue_type_set_name);
	stmt.bind(2,residue_type_name);
	cppdb::result res(basic::database::safely_read_from_database(stmt));
	while(res.next())
	{
		std::string property;
		res >> property;
		res_type.add_property(property);
	}
}



void
ResidueDatabaseIO::report_residue_type_variant(
	std::string const & residue_type_set_name,
	ResidueType const & res_type,
	utility::sql_database::sessionOP db_session
)  {
	foreach(std::string const & variant_type, res_type.variant_types()){

		cppdb::statement stmt = (*db_session)
			<< "INSERT INTO residue_type_variant_type VALUES (?,?,?);"
			<< residue_type_set_name
			<< res_type.name()
			<< variant_type;
		basic::database::safely_write_to_database(stmt);

	}
}

void
ResidueDatabaseIO::read_residue_type_variant(
	std::string const & residue_type_set_name,
	std::string const & residue_type_name,
	core::chemical::ResidueType & res_type,
	utility::sql_database::sessionOP db_session)
{
	std::string residue_variant_statement_string =
		"SELECT\n"
		"	variant_type\n"
		"FROM residue_type_variant_type\n"
		"WHERE residue_type_set_name = ? AND residue_type_name = ?;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(residue_variant_statement_string,db_session));
	stmt.bind(1,residue_type_set_name);
	stmt.bind(2,residue_type_name);
	cppdb::result res(basic::database::safely_read_from_database(stmt));
	while(res.next())
	{
		std::string variant_type;
		res >> variant_type;
		res_type.add_variant_type(variant_type);

	}
}

void
ResidueDatabaseIO::report_residue_type_icoor(
	std::string const & residue_type_set_name,
	core::chemical::ResidueType const & res_type,
	utility::sql_database::sessionOP db_session)
{
	std::string statement_string = "INSERT INTO residue_type_icoor VALUES (?,?,?,?,?,?,?,?,?,?);";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));

	for(Size i=1; i <= res_type.natoms(); ++i){

		AtomICoor atom_icoor(res_type.icoor(i));
		stmt.bind(1,residue_type_set_name);
		stmt.bind(2,res_type.name());
		stmt.bind(3,atom_icoor.index());
		stmt.bind(4,res_type.atom_name(i));
		stmt.bind(5,atom_icoor.phi());
		stmt.bind(6,atom_icoor.theta());
		stmt.bind(7,atom_icoor.d());
		//We need to do this because the AtomICoor doesn't store atom names, just numbers, but it takes atom names as its constructor
		if(atom_icoor.stub_atom1().type() == ICoorAtomID::INTERNAL)
		{
			stmt.bind(8,res_type.atom_name(atom_icoor.stub_atom1().atomno()));
		}else if(atom_icoor.stub_atom1().type() == ICoorAtomID::CONNECT)
		{
			stmt.bind(8,"CONNECT");
		}else if(atom_icoor.stub_atom2().type() == ICoorAtomID::POLYMER_LOWER)
		{
			stmt.bind(8,"LOWER");
		}else
		{
			//should be POLYMER_UPPER
			stmt.bind(8,"UPPER");
		}

		if(atom_icoor.stub_atom2().type() == ICoorAtomID::INTERNAL)
		{
			stmt.bind(9,res_type.atom_name(atom_icoor.stub_atom2().atomno()));
		}else if(atom_icoor.stub_atom1().type() == ICoorAtomID::CONNECT)
		{
			stmt.bind(9,"CONNECT");
		}else if(atom_icoor.stub_atom2().type() == ICoorAtomID::POLYMER_LOWER)
		{
			stmt.bind(9,"LOWER");
		}else
		{
			//should be POLYMER_UPPER
			stmt.bind(9,"UPPER");
		}

		if(atom_icoor.stub_atom3().type() == ICoorAtomID::INTERNAL)
		{
			stmt.bind(10,res_type.atom_name(atom_icoor.stub_atom3().atomno()));
		}else if(atom_icoor.stub_atom1().type() == ICoorAtomID::CONNECT)
		{
			stmt.bind(10,"CONNECT");
		}else if(atom_icoor.stub_atom2().type() == ICoorAtomID::POLYMER_LOWER)
		{
			stmt.bind(10,"LOWER");
		}else
		{
			//should be POLYMER_UPPER
			stmt.bind(10,"UPPER");
		}
		basic::database::safely_write_to_database(stmt);

	}
}

void
ResidueDatabaseIO::read_residue_type_icoor(
	std::string const & residue_type_set_name,
	std::string const & residue_type_name,
	core::chemical::ResidueType & res_type,
	utility::sql_database::sessionOP db_session)
{

	//ORDER BY icoor_sequence exists because the icoor assignment code is order dependent
	std::string residue_variant_statement_string =
		"SELECT\n"
		"	child_atom,\n"
		"	phi,\n"
		"	theta,\n"
		"	distance,\n"
		"	parent_atom,\n"
		"	angle_atom,\n"
		"	torsion_atom\n"
		"FROM residue_type_icoor\n"
		"WHERE residue_type_set_name = ? AND residue_type_name = ?\n"
		"ORDER BY icoor_sequence;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(residue_variant_statement_string,db_session));
	stmt.bind(1,residue_type_set_name);
	stmt.bind(2,residue_type_name);

	std::map< std::string, Vector > rsd_xyz;

	cppdb::result res(basic::database::safely_read_from_database(stmt));
	while(res.next())
	{
		std::string child_atom_name, parent_atom_name,angle_atom_name,torsion_atom_name;
		core::Real phi, theta, distance;
		res >> child_atom_name >> phi >> theta >> distance >> parent_atom_name >>angle_atom_name >>torsion_atom_name;

		//This code mostly came from core/chemical/residue_io.  It should probably be extracted to a util function
		if(child_atom_name == parent_atom_name)
		{
			assert( rsd_xyz.empty() ); // first atom
			rsd_xyz[child_atom_name] = Vector(0.0);
		}else if( child_atom_name == angle_atom_name)
		{
			assert( rsd_xyz.size() == 1 && rsd_xyz.count( parent_atom_name ) ); // second atom
			rsd_xyz[ child_atom_name ] = Vector( distance, 0.0, 0.0 );
		}else
		{
			Vector torsion_xyz;
			if ( child_atom_name == torsion_atom_name ) {
				assert( rsd_xyz.size() == 2 );
				assert( rsd_xyz.count( parent_atom_name ) );
				assert( rsd_xyz.count( angle_atom_name ) ); // third atom
				torsion_xyz = Vector( 1.0, 1.0, 0.0 );
			} else {
				assert( rsd_xyz.count( parent_atom_name ) && rsd_xyz.count( angle_atom_name ) && rsd_xyz.count( torsion_atom_name ) );
				torsion_xyz = rsd_xyz[ torsion_atom_name ];
			}
			kinematics::Stub const stub( rsd_xyz[ parent_atom_name ], rsd_xyz[ angle_atom_name ], torsion_xyz );
			rsd_xyz[ child_atom_name ] = stub.spherical( phi, theta, distance );
		}

		// set atom_base
		if ( child_atom_name != "UPPER" && child_atom_name != "LOWER" && child_atom_name.substr(0,4) != "CONN" ) {
			// atom base only valid for genuine atoms of this residue
			if ( child_atom_name == parent_atom_name ) {
				// root of the tree
				if ( res_type.natoms() == 1 ) {
					res_type.set_atom_base( child_atom_name, child_atom_name ); // 1st child of root atom
				} else {
					res_type.set_atom_base( child_atom_name, angle_atom_name ); // 1st child of root atom
				}
			} else {
				res_type.set_atom_base( child_atom_name, parent_atom_name );
			}
		}

		// set icoor
		res_type.set_icoor( child_atom_name, phi, theta, distance, parent_atom_name, angle_atom_name, torsion_atom_name );
	}

	for ( Size i=1; i<= res_type.natoms(); ++i ) {
		std::string name( res_type.atom_name(i) );
		assert( rsd_xyz.count( name ) );
		res_type.set_xyz( name, rsd_xyz[ name ] );
	}


}




}
}
