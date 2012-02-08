// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/chemical/AtomTypeDatabaseIO.cc
/// @brief Write chemical information associated with atom types to a database
/// @author Matthew O'Meara

// Project Headers
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomTypeDatabaseIO.hh>
#include <core/chemical/types.hh>
#include <basic/database/sql_utils.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>

#include <utility/vector1.hh>

// Boost Headers
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

namespace core {
namespace chemical {

using std::string;
using core::Real;

using utility::vector1;
using core::chemical::AtomType;
using utility::sql_database::sessionOP;
using basic::database::write_schema_to_database;
using basic::database::safely_prepare_statement;
using basic::database::safely_read_from_database;
using basic::database::safely_write_to_database;

using cppdb::statement;
using cppdb::result;

AtomTypeDatabaseIO::AtomTypeDatabaseIO() {}

AtomTypeDatabaseIO::~AtomTypeDatabaseIO() {}


//////////////////////////
//                      //
//    Public Interface  //
//                      //
//////////////////////////

std::string AtomTypeDatabaseIO::schema() const
{
	// NOTE: To support building feature databases dentified by
	// their names rather then assigning them a unique id.

	string db_mode(basic::options::option[basic::options::OptionKeys::inout::database_mode]);
	if(db_mode == "sqlite3") {
		return
			"CREATE TABLE IF NOT EXISTS atom_types (\n"
			"	atom_type_set_name TEXT,\n"
			"	name TEXT,\n"
			"	element TEXT,\n"
			"	lennard_jones_radius REAL,\n"
			"	lennard_jones_well_depth REAL,\n"
			"	lazaridis_karplus_lambda REAL,\n"
			"	lazaridis_karplus_degrees_of_freedom REAL,\n"
			"	lazaridis_karplus_volume REAL,\n"
			"	PRIMARY KEY(atom_type_set_name, name));\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS atom_type_property_values (\n"
			"	property TEXT,\n"
			"	PRIMARY KEY(property));\n"
			"\n"
			"INSERT OR IGNORE INTO atom_type_property_values VALUES ('ACCEPTOR');\n"
			"INSERT OR IGNORE INTO atom_type_property_values VALUES ('DONOR');\n"
			"INSERT OR IGNORE INTO atom_type_property_values VALUES ('POLAR_HYDROGEN');\n"
			"INSERT OR IGNORE INTO atom_type_property_values VALUES ('AROMATIC');\n"
			"INSERT OR IGNORE INTO atom_type_property_values VALUES ('H2O');\n"
			"INSERT OR IGNORE INTO atom_type_property_values VALUES ('ORBITALS');\n"
			"INSERT OR IGNORE INTO atom_type_property_values VALUES ('VIRTUAL');\n"
			"INSERT OR IGNORE INTO atom_type_property_values VALUES ('SP2_HYBRID');\n"
			"INSERT OR IGNORE INTO atom_type_property_values VALUES ('SP3_HYBRID');\n"
			"INSERT OR IGNORE INTO atom_type_property_values VALUES ('RING_HYBRID');\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS atom_type_properties (\n"
			"	atom_type_set_name TEXT,\n"
			"	name TEXT,\n"
			"	property TEXT,\n"
			"	FOREIGN KEY(atom_type_set_name, name)\n"
			"		REFERENCES atom_types (atom_type_set_name, name)\n"
			"		DEFERRABLE INITIALLY DEFERRED,\n"
			"	FOREIGN KEY(property)\n"
			"		REFERENCES atom_type_property_values (property)\n"
			"		DEFERRABLE INITIALLY DEFERRED,\n"
			"	PRIMARY KEY(atom_type_set_name, name, property));\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS atom_type_extra_parameters (\n"
			"	atom_type_set_name TEXT,\n"
			"	name TEXT,\n"
			"	parameter TEXT,\n"
			"	value REAL,\n"
			"	FOREIGN KEY(atom_type_set_name, name)\n"
			"		REFERENCES atom_types (atom_type_set_name, name)\n"
			"		DEFERRABLE INITIALLY DEFERRED,\n"
			"	PRIMARY KEY(atom_type_set_name, name, parameter));";
	} else if(db_mode == "mysql") {
		return
			"CREATE TABLE IF NOT EXISTS atom_types (\n"
			"	atom_type_set_name VARCHAR(64),\n"
			"	name VARCHAR(4),\n"
			"	element VARCHAR(1),\n"
			"	lennard_jones_radius REAL,\n"
			"	lennard_jones_well_depth REAL,\n"
			"	lazaridis_karplus_lambda REAL,\n"
			"	lazaridis_karplus_degrees_of_freedom REAL,\n"
			"	lazaridis_karplus_volume REAL,\n"
			"	PRIMARY KEY(atom_type_set_name, name));\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS atom_type_proproperty_values (\n"
			"	property VARCHAR(32),\n"
			"	PRIMARY KEY(property));\n"
			"\n"
			"INSERT OR IGNORE INTO atom_type_property_values VALUES ( 'ACCEPTOR' );\n"
			"INSERT OR IGNORE INTO atom_type_property_values VALUES ( 'DONOR' );\n"
			"INSERT OR IGNORE INTO atom_type_property_values VALUES ( 'POLAR_HYDROGEN' );\n"
			"INSERT OR IGNORE INTO atom_type_property_values VALUES ( 'AROMATIC' );\n"
			"INSERT OR IGNORE INTO atom_type_property_values VALUES ( 'H2O' );\n"
			"INSERT OR IGNORE INTO atom_type_property_values VALUES ( 'ORBITALS' );\n"
			"INSERT OR IGNORE INTO atom_type_property_values VALUES ( 'VIRTUAL' );\n"
			"INSERT OR IGNORE INTO atom_type_property_values VALUES ( 'SP2_HYBRID' );\n"
			"INSERT OR IGNORE INTO atom_type_property_values VALUES ( 'SP3_HYBRID' );\n"
			"INSERT OR IGNORE INTO atom_type_property_values VALUES ( 'RING_HYBRID' );\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS atom_type_properties (\n"
			"	atom_type_set_name VARCHAR(64),\n"
			"	name (4),\n"
			"	property VARCHAR(32),\n"
			"	FOREIGN KEY(atom_type_set_name, name)\n"
			"		REFERENCES atom_types (atom_type_set_name, name),\n"
			"	FOREIGN KEY(property)\n"
			"		REFERENCES atom_type_property_values (property),\n"
			"	PRIMARY KEY(atom_type_set_name, name, property))\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS atom_type_extra_parameters (\n"
			"	atom_type_set_name VARCHAR(64),\n"
			"	name (4),\n"
			"	parameter TEXT,\n"
			"	value REAL,\n"
			"	FOREIGN KEY(atom_type_set_name, name)\n"
			"		REFERENCES atom_types (atom_type_set_name, name)\n"
			"	PRIMARY KEY(atom_type_set_name, name, parameter));";
	} else {
		return "";
	}

}

void
AtomTypeDatabaseIO::initialize(
	sessionOP db_session
) const {
	write_schema_to_database(schema(),db_session);
}

void
AtomTypeDatabaseIO::write_atom_type_set_to_database(
	AtomTypeSet const & atom_type_set,
	sessionOP db_session
) const {
	string const & atom_type_set_name(atom_type_set.name());

	string stmt_string =
		"SELECT * FROM atom_types WHERE atom_type_set_name = ?";
	statement stmt(safely_prepare_statement(stmt_string, db_session));
	stmt.bind(1, atom_type_set_name);
	result res(safely_read_from_database(stmt));
	if(res.next()) return;


	for(Size atom_index=1; atom_index <= atom_type_set.n_atomtypes(); ++atom_index){
		AtomType const & atom_type(atom_type_set[atom_index]);
		write_atom_type_table(
			atom_type_set_name, atom_type, db_session);
		write_atom_type_properties_table(
			atom_type_set_name, atom_type, db_session);
		write_atom_type_extra_parameters_table(
			atom_type_set, atom_type, db_session);
	}
}

vector1<string>
AtomTypeDatabaseIO::get_all_atom_types_in_database(
	sessionOP db_session
) const {

	vector1<string> atom_names;

	string stmt_string = "SELECT name FROM atom_types;";
	statement stmt(safely_prepare_statement(stmt_string, db_session));
	result res(safely_read_from_database(stmt));
	while(res.next())
	{
		string name;
		res >> name;
		atom_names.push_back(name);
	}
	return atom_names;
}


///////////////////////
//                   //
//  Helper Methods   //
//                   //
///////////////////////


void
AtomTypeDatabaseIO::write_atom_type_table(
	string const & atom_type_set_name,
	AtomType const & atom_type,
	sessionOP db_session
) const {

	string stmt_string = "INSERT INTO atom_types VALUES (?,?,?,?,?,?,?,?);";
	statement stmt(safely_prepare_statement(stmt_string, db_session));

	stmt.bind(1,atom_type_set_name);
	stmt.bind(2,atom_type.name());
	stmt.bind(3,atom_type.element());
	stmt.bind(4,atom_type.lj_radius());
	stmt.bind(5,atom_type.lj_wdepth());
	stmt.bind(6,atom_type.lk_lambda());
	stmt.bind(7,atom_type.lk_dgfree());
	stmt.bind(8,atom_type.lk_volume());
	safely_write_to_database(stmt);

}


void
AtomTypeDatabaseIO::write_atom_type_properties_table(
	string const & atom_type_set_name,
	AtomType const & atom_type,
	sessionOP db_session
) const {

	string statement_string = "INSERT INTO atom_type_properties VALUES (?,?,?);";
	statement stmt(safely_prepare_statement(statement_string, db_session));

	vector1<string> properties(atom_type.get_all_properties());
	for(vector1<string>::const_iterator
			property_iter = properties.begin(),	property_end = properties.end();
			property_iter != property_end; ++property_iter){
		stmt.bind(1, atom_type_set_name);
		stmt.bind(2, atom_type.name());
		stmt.bind(3, *property_iter);
		safely_write_to_database(stmt);
	}
}


void
AtomTypeDatabaseIO::write_atom_type_extra_parameters_table(
	AtomTypeSet const & atom_type_set,
	AtomType const & atom_type,
	sessionOP db_session
) const {

	string stmt_string = "INSERT INTO atom_type_extra_parameters VALUES (?,?,?,?);";
	statement stmt(safely_prepare_statement(stmt_string, db_session));

	for(std::map<std::string, int>::const_iterator
		extra_parameter_index_iter =
			atom_type_set.extra_parameter_indices().begin(),
		extra_parameter_index_iter_end =
			atom_type_set.extra_parameter_indices().end();
		extra_parameter_index_iter != extra_parameter_index_iter_end;
		++extra_parameter_index_iter){

		string const extra_parameter_name(extra_parameter_index_iter->first);
		Size const extra_parameter_index(extra_parameter_index_iter->second);
		stmt.bind(1, atom_type_set.name());
		stmt.bind(2, atom_type.name());
		stmt.bind(3, extra_parameter_name);
		stmt.bind(4, atom_type.extra_parameter(extra_parameter_index));
		safely_write_to_database(stmt);
	}
}


} // namespace
} // namespace
