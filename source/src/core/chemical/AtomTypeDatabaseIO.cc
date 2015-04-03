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
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>


#include <utility/vector1.hh>

// Boost Headers
#include <boost/foreach.hpp>
#include <boost/assign/list_of.hpp>

#include <vector>

namespace core {
namespace chemical {

using std::string;
using core::Real;

using utility::vector1;
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


void
AtomTypeDatabaseIO::write_schema_to_db(
	sessionOP db_session
) const {
	write_atom_types_table_schema(db_session);
	write_atom_type_property_values_table_schema(db_session);
	write_atom_type_properties_table_schema(db_session);
	write_atom_type_extra_parameters_table_schema(db_session);
}

void
AtomTypeDatabaseIO::write_atom_types_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column atom_type_set_name("atom_type_set_name", DbDataTypeOP( new DbText(64) ));
	Column name("name", DbDataTypeOP( new DbText(32) ));
	Column element("element", DbDataTypeOP( new DbText(2) ));
	Column lennard_jones_radius("lennard_jones_radius", DbDataTypeOP( new DbReal() ));
	Column lennard_jones_well_depth("lennard_jones_well_depth", DbDataTypeOP( new DbReal() ));
	Column lazaridis_karplus_lambda("lazaridis_karplus_lambda", DbDataTypeOP( new DbReal() ));
	Column lazaridis_karplus_degrees_of_freedom("lazaridis_karplus_degrees_of_freedom", DbDataTypeOP( new DbReal() ));
	Column lazaridis_karplus_volume("lazaridis_karplus_volume", DbDataTypeOP( new DbReal() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(atom_type_set_name);
	primary_key_columns.push_back(name);
	PrimaryKey primary_key(primary_key_columns);

	Schema table("atom_types", primary_key);
	table.add_column(element);
	table.add_column(lennard_jones_radius);
	table.add_column(lennard_jones_well_depth);
	table.add_column(lazaridis_karplus_lambda);
	table.add_column(lazaridis_karplus_degrees_of_freedom);
	table.add_column(lazaridis_karplus_volume);

	table.write(db_session);
}

void
AtomTypeDatabaseIO::write_atom_type_property_values_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;
	using namespace basic::database;
	using namespace boost::assign;

	Column property("property", DbDataTypeOP( new DbText(32) ));

	Columns primary_key_columns;
	primary_key_columns.push_back(property);
	PrimaryKey primary_key(primary_key_columns);

	Schema table("atom_type_property_values", primary_key);

	table.write(db_session);

	// insert values
	string const table_name("atom_type_property_values");
	std::vector<string> column_names;
	column_names.push_back("property");
	insert_or_ignore(table_name, column_names, list_of("ACCEPTOR"), db_session);
	insert_or_ignore(table_name, column_names, list_of("ACCEPTOR"), db_session);
	insert_or_ignore(table_name, column_names, list_of("DONOR"), db_session);
	insert_or_ignore(table_name, column_names, list_of("POLAR_HYDROGEN"), db_session);
	insert_or_ignore(table_name, column_names, list_of("AROMATIC"), db_session);
	insert_or_ignore(table_name, column_names, list_of("H2O"), db_session);
	insert_or_ignore(table_name, column_names, list_of("ORBITALS"), db_session);
	insert_or_ignore(table_name, column_names, list_of("VIRTUAL"), db_session);
	insert_or_ignore(table_name, column_names, list_of("SP2_HYBRID"), db_session);
	insert_or_ignore(table_name, column_names, list_of("SP3_HYBRID"), db_session);
	insert_or_ignore(table_name, column_names, list_of("RING_HYBRID"), db_session);

}

void
AtomTypeDatabaseIO::write_atom_type_properties_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column atom_type_set_name("atom_type_set_name", DbDataTypeOP( new DbText(64) ));
	Column name("name", DbDataTypeOP( new DbText(4) ));
	Column property("property", DbDataTypeOP( new DbText(32) ));

	Columns primary_key_columns;
	primary_key_columns.push_back(atom_type_set_name);
	primary_key_columns.push_back(name);
	primary_key_columns.push_back(property);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns1;
	foreign_key_columns1.push_back(atom_type_set_name);
	foreign_key_columns1.push_back(name);
	vector1< std::string > reference_columns1;
	reference_columns1.push_back("atom_type_set_name");
	reference_columns1.push_back("name");
	ForeignKey foreign_key1(foreign_key_columns1, "atom_types", reference_columns1, true);

	Columns foreign_key_columns2;
	foreign_key_columns2.push_back(property);
	vector1< std::string > reference_columns2;
	reference_columns2.push_back("property");
	ForeignKey foreign_key2(foreign_key_columns2, "atom_type_property_values", reference_columns2, true);


	Schema table("atom_type_properties", primary_key);
	table.add_foreign_key(foreign_key1);
	table.add_foreign_key(foreign_key2);

	table.write(db_session);
}


void
AtomTypeDatabaseIO::write_atom_type_extra_parameters_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column atom_type_set_name("atom_type_set_name", DbDataTypeOP( new DbText(64) ));
	Column name("name", DbDataTypeOP( new DbText(32) ));
	Column parameter("parameter", DbDataTypeOP( new DbText(32) ));
	Column value("value", DbDataTypeOP( new DbReal() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(atom_type_set_name);
	primary_key_columns.push_back(name);
	primary_key_columns.push_back(parameter);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(atom_type_set_name);
	foreign_key_columns.push_back(name);
	vector1< std::string > reference_columns;
	reference_columns.push_back("atom_type_set_name");
	reference_columns.push_back("name");
	ForeignKey foreign_key(foreign_key_columns, "atom_types", reference_columns, true);

	Schema table("atom_type_extra_parameters", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(value);

	table.write(db_session);
}


void
AtomTypeDatabaseIO::initialize(
	sessionOP db_session
) const {
	write_schema_to_db(db_session);
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

	string stmt_string = "INSERT INTO atom_types (atom_type_set_name, name, element, lennard_jones_radius, lennard_jones_well_depth, lazaridis_karplus_lambda, lazaridis_karplus_degrees_of_freedom, lazaridis_karplus_volume) VALUES (?,?,?,?,?,?,?,?);";
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

	string statement_string = "INSERT INTO atom_type_properties (atom_type_set_name, name, property) VALUES (?,?,?);";
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

	string const atom_type_set_name = atom_type_set.name();

	string stmt_string = "INSERT INTO atom_type_extra_parameters (atom_type_set_name, name, parameter, value) VALUES (?,?,?,?);";
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
		stmt.bind(1, atom_type_set_name);
		stmt.bind(2, atom_type.name());
		stmt.bind(3, extra_parameter_name);
		stmt.bind(4, atom_type.extra_parameter(extra_parameter_index));
		safely_write_to_database(stmt);
	}
}


} // namespace
} // namespace
