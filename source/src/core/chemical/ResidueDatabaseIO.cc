// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/chemical/ResidueDatabaseIO.cc
/// @author Sam DeLuca
/// @author Matt O'Meara

// Project Headers
#include <core/chemical/Atom.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/ResidueDatabaseIO.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/icoor_support.hh>
#include <core/chemical/rotamers/PDBRotamerLibrarySpecification.hh>
#include <core/types.hh>
#include <basic/database/sql_utils.hh>

#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/Index.hh>
#include <basic/database/schema_generator/DbDataType.hh>

#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>

// Boost Headers
#include <boost/foreach.hpp>

namespace core {
namespace chemical {

using utility::sql_database::sessionOP;
using utility::vector1;

ResidueDatabaseIO::ResidueDatabaseIO() :version_(0.2)
{

}

ResidueDatabaseIO::~ResidueDatabaseIO() = default;

void
ResidueDatabaseIO::write_schema_to_db(
	sessionOP db_session
) const {
	write_residue_type_table_schema(db_session);
	write_residue_type_atom_table_schema(db_session);
	write_residue_type_bond_table_schema(db_session);
	write_residue_type_cut_bond_table_schema(db_session);
	write_residue_type_chi_table_schema(db_session);
	write_residue_type_chi_rotamer_table_schema(db_session);
	write_residue_type_proton_chi_table_schema(db_session);
	write_residue_type_property_table_schema(db_session);
	write_residue_type_variant_type_table_schema(db_session);
	write_residue_type_icoor_table_schema(db_session);
}

void
ResidueDatabaseIO::write_residue_type_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column residue_type_set_name("residue_type_set_name", DbDataTypeOP( new DbText(255) ));
	Column version("version", DbDataTypeOP( new DbText() ));
	Column name("name", DbDataTypeOP( new DbText(255) ));
	Column name3("name3", DbDataTypeOP( new DbText(3) ));
	Column name1("name1", DbDataTypeOP( new DbText(2) )); // TODO Fix ccpdb to allow putting in just a single character rather then a one character string.
	Column aa("aa", DbDataTypeOP( new DbInteger() ));
	Column lower_connect("lower_connect", DbDataTypeOP( new DbInteger() ));
	Column upper_connect("upper_connect", DbDataTypeOP( new DbInteger() ));
	Column nbr_atom("nbr_atom", DbDataTypeOP( new DbInteger() ));
	Column nbr_radius("nbr_radius", DbDataTypeOP( new DbReal() ));
	Column rotamer_library("rotamer_library", DbDataTypeOP( new DbText() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(residue_type_set_name);
	primary_key_columns.push_back(name);
	PrimaryKey primary_key(primary_key_columns);

	Schema table("residue_type", primary_key);
	table.add_column(version);
	table.add_column(name3);
	table.add_column(name1);
	table.add_column(aa);
	table.add_column(lower_connect);
	table.add_column(upper_connect);
	table.add_column(nbr_atom);
	table.add_column(nbr_radius);
	table.add_column(rotamer_library);

	table.write(db_session);
}

void
ResidueDatabaseIO::write_residue_type_atom_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column residue_type_set_name("residue_type_set_name", DbDataTypeOP( new DbText(255) ));
	Column residue_type_name("residue_type_name", DbDataTypeOP( new DbText(255) ));
	Column atom_index("atom_index", DbDataTypeOP( new DbInteger() ));
	Column atom_name("atom_name", DbDataTypeOP( new DbText() ));
	Column atom_type_name("atom_type_name", DbDataTypeOP( new DbText() ));
	Column mm_atom_type_name("mm_atom_type_name", DbDataTypeOP( new DbText() ));
	Column charge("charge", DbDataTypeOP( new DbReal() ));
	Column is_backbone("is_backbone", DbDataTypeOP( new DbInteger() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(residue_type_set_name);
	primary_key_columns.push_back(residue_type_name);
	primary_key_columns.push_back(atom_index);
	PrimaryKey primary_key(primary_key_columns);

	Columns index1_columns;
	index1_columns.push_back(residue_type_set_name);
	index1_columns.push_back(residue_type_name);
	index1_columns.push_back(atom_name);
	Index index1(index1_columns);


	Columns foreign_key_columns;
	foreign_key_columns.push_back(residue_type_set_name);
	foreign_key_columns.push_back(residue_type_name);
	vector1< std::string > reference_columns;
	reference_columns.push_back("residue_type_set_name");
	reference_columns.push_back("name");
	ForeignKey foreign_key(foreign_key_columns, "residue_type", reference_columns, true);

	Schema table("residue_type_atom", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(atom_name);
	table.add_column(atom_type_name);
	table.add_column(mm_atom_type_name);
	table.add_column(charge);
	table.add_column(is_backbone);
	table.add_index(index1);

	table.write(db_session);
}

void
ResidueDatabaseIO::write_residue_type_bond_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column residue_type_set_name("residue_type_set_name", DbDataTypeOP( new DbText(255) ));
	Column residue_type_name("residue_type_name", DbDataTypeOP( new DbText(255) ));
	Column atom1("atom1", DbDataTypeOP( new DbInteger() ));
	Column atom2("atom2", DbDataTypeOP( new DbInteger() ));
	Column bond_type("bond_type", DbDataTypeOP( new DbInteger() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(residue_type_set_name);
	primary_key_columns.push_back(residue_type_name);
	primary_key_columns.push_back(atom1);
	primary_key_columns.push_back(atom2);
	PrimaryKey primary_key(primary_key_columns);


	Columns foreign_key_columns;
	foreign_key_columns.push_back(residue_type_set_name);
	foreign_key_columns.push_back(residue_type_name);
	vector1< std::string > reference_columns;
	reference_columns.push_back("residue_type_set_name");
	reference_columns.push_back("name");
	ForeignKey foreign_key(foreign_key_columns, "residue_type", reference_columns, true);

	Schema table("residue_type_bond", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(bond_type);

	table.write(db_session);
}


void
ResidueDatabaseIO::write_residue_type_cut_bond_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column residue_type_set_name("residue_type_set_name", DbDataTypeOP( new DbText(255) ));
	Column residue_type_name("residue_type_name", DbDataTypeOP( new DbText(255) ));
	Column atom1("atom1", DbDataTypeOP( new DbInteger() ));
	Column atom2("atom2", DbDataTypeOP( new DbInteger() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(residue_type_set_name);
	primary_key_columns.push_back(residue_type_name);
	primary_key_columns.push_back(atom1);
	primary_key_columns.push_back(atom2);
	PrimaryKey primary_key(primary_key_columns);


	Columns foreign_key_columns;
	foreign_key_columns.push_back(residue_type_set_name);
	foreign_key_columns.push_back(residue_type_name);
	vector1< std::string > reference_columns;
	reference_columns.push_back("residue_type_set_name");
	reference_columns.push_back("name");
	ForeignKey foreign_key(foreign_key_columns, "residue_type", reference_columns, true);

	Schema table("residue_type_cut_bond", primary_key);
	table.add_foreign_key(foreign_key);

	table.write(db_session);
}

void
ResidueDatabaseIO::write_residue_type_chi_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column residue_type_set_name("residue_type_set_name", DbDataTypeOP( new DbText(255) ));
	Column residue_type_name("residue_type_name", DbDataTypeOP( new DbText(255) ));
	Column chino("chino", DbDataTypeOP( new DbInteger() ));
	Column atom1("atom1", DbDataTypeOP( new DbText() ));
	Column atom2("atom2", DbDataTypeOP( new DbText() ));
	Column atom3("atom3", DbDataTypeOP( new DbText() ));
	Column atom4("atom4", DbDataTypeOP( new DbText() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(residue_type_set_name);
	primary_key_columns.push_back(residue_type_name);

	// perhaps it would be better to have chino here instead of the
	// atoms themselves?
	primary_key_columns.push_back(atom1);
	primary_key_columns.push_back(atom2);
	primary_key_columns.push_back(atom3);
	primary_key_columns.push_back(atom4);
	PrimaryKey primary_key(primary_key_columns);


	Columns foreign_key_columns;
	foreign_key_columns.push_back(residue_type_set_name);
	foreign_key_columns.push_back(residue_type_name);
	vector1< std::string > reference_columns;
	reference_columns.push_back("residue_type_set_name");
	reference_columns.push_back("name");
	ForeignKey foreign_key(foreign_key_columns, "residue_type", reference_columns, true);

	Schema table("residue_type_chi", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(chino);

	table.write(db_session);
}

void
ResidueDatabaseIO::write_residue_type_chi_rotamer_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column residue_type_set_name("residue_type_set_name", DbDataTypeOP( new DbText(255) ));
	Column residue_type_name("residue_type_name", DbDataTypeOP( new DbText(255) ));
	Column chino("chino", DbDataTypeOP( new DbInteger() ));
	Column mean("mean", DbDataTypeOP( new DbReal() ));
	Column sdev("sdev", DbDataTypeOP( new DbReal() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(residue_type_set_name);
	primary_key_columns.push_back(residue_type_name);
	primary_key_columns.push_back(chino);
	primary_key_columns.push_back(mean);
	primary_key_columns.push_back(sdev);
	PrimaryKey primary_key(primary_key_columns);


	Columns foreign_key_columns;
	foreign_key_columns.push_back(residue_type_set_name);
	foreign_key_columns.push_back(residue_type_name);
	vector1< std::string > reference_columns;
	reference_columns.push_back("residue_type_set_name");
	reference_columns.push_back("name");
	ForeignKey foreign_key(foreign_key_columns, "residue_type", reference_columns, true);

	Schema table("residue_type_chi_rotamer", primary_key);
	table.add_foreign_key(foreign_key);

	table.write(db_session);
}

void
ResidueDatabaseIO::write_residue_type_proton_chi_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column residue_type_set_name("residue_type_set_name", DbDataTypeOP( new DbText(255) ));
	Column residue_type_name("residue_type_name", DbDataTypeOP( new DbText(255) ));
	Column chino("chino", DbDataTypeOP( new DbInteger() ));
	Column sample("sample", DbDataTypeOP( new DbReal() ));
	Column is_extra("is_extra", DbDataTypeOP( new DbInteger() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(residue_type_set_name);
	primary_key_columns.push_back(residue_type_name);
	primary_key_columns.push_back(chino);
	primary_key_columns.push_back(sample);
	PrimaryKey primary_key(primary_key_columns);


	Columns foreign_key_columns;
	foreign_key_columns.push_back(residue_type_set_name);
	foreign_key_columns.push_back(residue_type_name);
	vector1< std::string > reference_columns;
	reference_columns.push_back("residue_type_set_name");
	reference_columns.push_back("name");
	ForeignKey foreign_key(foreign_key_columns, "residue_type", reference_columns, true);

	Schema table("residue_type_proton_chi", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(is_extra);

	table.write(db_session);
}

void
ResidueDatabaseIO::write_residue_type_property_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column residue_type_set_name("residue_type_set_name", DbDataTypeOP( new DbText(255) ));
	Column residue_type_name("residue_type_name", DbDataTypeOP( new DbText(255) ));
	Column property("property", DbDataTypeOP( new DbText() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(residue_type_set_name);
	primary_key_columns.push_back(residue_type_name);
	primary_key_columns.push_back(property);
	PrimaryKey primary_key(primary_key_columns);


	Columns foreign_key_columns;
	foreign_key_columns.push_back(residue_type_set_name);
	foreign_key_columns.push_back(residue_type_name);
	vector1< std::string > reference_columns;
	reference_columns.push_back("residue_type_set_name");
	reference_columns.push_back("name");
	ForeignKey foreign_key(foreign_key_columns, "residue_type", reference_columns, true);

	Schema table("residue_type_property", primary_key);
	table.add_foreign_key(foreign_key);

	table.write(db_session);
}

void
ResidueDatabaseIO::write_residue_type_variant_type_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column residue_type_set_name("residue_type_set_name", DbDataTypeOP( new DbText(255) ));
	Column residue_type_name("residue_type_name", DbDataTypeOP( new DbText(255) ));
	Column variant_type("variant_type", DbDataTypeOP( new DbText() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(residue_type_set_name);
	primary_key_columns.push_back(residue_type_name);
	primary_key_columns.push_back(variant_type);
	PrimaryKey primary_key(primary_key_columns);


	Columns foreign_key_columns;
	foreign_key_columns.push_back(residue_type_set_name);
	foreign_key_columns.push_back(residue_type_name);
	vector1< std::string > reference_columns;
	reference_columns.push_back("residue_type_set_name");
	reference_columns.push_back("name");
	ForeignKey foreign_key(foreign_key_columns, "residue_type", reference_columns, true);

	Schema table("residue_type_variant_type", primary_key);
	table.add_foreign_key(foreign_key);

	table.write(db_session);
}

void
ResidueDatabaseIO::write_residue_type_icoor_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column residue_type_set_name("residue_type_set_name", DbDataTypeOP( new DbText(255) ));
	Column residue_type_name("residue_type_name", DbDataTypeOP( new DbText(255) ));
	Column icoor_sequence("icoor_sequence", DbDataTypeOP( new DbInteger() ));
	Column child_atom("child_atom", DbDataTypeOP( new DbText() ));
	Column phi("phi", DbDataTypeOP( new DbReal() ));
	Column theta("theta", DbDataTypeOP( new DbReal() ));
	Column distance("distance", DbDataTypeOP( new DbReal() ));
	Column parent_atom("parent_atom", DbDataTypeOP( new DbText() ));
	Column angle_atom("angle_atom", DbDataTypeOP( new DbText() ));
	Column torsion_atom("torsion_atom", DbDataTypeOP( new DbText() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(residue_type_set_name);
	primary_key_columns.push_back(residue_type_name);
	primary_key_columns.push_back(child_atom);
	PrimaryKey primary_key(primary_key_columns);


	Columns foreign_key_columns;
	foreign_key_columns.push_back(residue_type_set_name);
	foreign_key_columns.push_back(residue_type_name);
	vector1< std::string > reference_columns;
	reference_columns.push_back("residue_type_set_name");
	reference_columns.push_back("name");
	ForeignKey foreign_key(foreign_key_columns, "residue_type", reference_columns, true);

	Schema table("residue_type_icoor", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(icoor_sequence);
	table.add_column(phi);
	table.add_column(theta);
	table.add_column(distance);
	table.add_column(parent_atom);
	table.add_column(angle_atom);
	table.add_column(torsion_atom);

	table.write(db_session);
}

void ResidueDatabaseIO::initialize(utility::sql_database::sessionOP db_session)
{
	write_schema_to_db(db_session);
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

	if ( res.next() ) return;

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
	chemical::AtomTypeSetCOP atom_types,
	chemical::ElementSetCOP elements,
	chemical::MMAtomTypeSetCOP mm_atom_types,
	chemical::orbitals::OrbitalTypeSetCOP orbital_atom_types,
	std::string const & residue_type_set_name,
	std::string const & residue_type_name,
	utility::sql_database::sessionOP db_session
)
{

	core::chemical::ResidueTypeOP res_type( new core::chemical::ResidueType(atom_types,elements,mm_atom_types,orbital_atom_types) );


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

	if ( atom_name_id_cache_.find(atom_id) != atom_name_id_cache_.end() ) {
		return atom_name_id_cache_[atom_id];
	}

	std::string atom_query =
		"SELECT\n"
		"\tatom_index,\n"
		"\tatom_name\n"
		"FROM residue_type_atom WHERE residue_type_name=?;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(atom_query,db_session));
	stmt.bind(1,residue_name);
	cppdb::result res(basic::database::safely_read_from_database(stmt));

	while ( res.next() )
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

	// This needs to be in a local variable as bind takes a reference, and the reference needs to be valid until the safely_write_to_database()
	std::stringstream name1_stream;
	name1_stream << res_type.name1();
	std::string name1( name1_stream.str() );

	int lower_terminus(-1), upper_terminus(-1);
	if ( res_type.is_polymer() ) {
		if ( !res_type.is_lower_terminus() ) lower_terminus = res_type.lower_connect_atom();
		if ( !res_type.is_upper_terminus() ) upper_terminus = res_type.upper_connect_atom();
	}

	std::string statement_string = "INSERT INTO residue_type (residue_type_set_name, name, version, name3, name1, aa, lower_connect, upper_connect, nbr_atom, nbr_radius, rotamer_library) VALUES (?,?,?,?,?,?,?,?,?,?,?);";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));

	stmt.bind(1,residue_type_set_name);
	stmt.bind(2,res_type.name());
	stmt.bind(3,version_);
	stmt.bind(4,res_type.name3());
	stmt.bind(5,name1);
	stmt.bind(6,res_type.aa());
	stmt.bind(7,lower_terminus);
	stmt.bind(8,upper_terminus);
	stmt.bind(9,res_type.nbr_atom());
	stmt.bind(10,res_type.nbr_radius());
	std::string libname("");
	rotamers::PDBRotamerLibrarySpecificationCOP pdb_rotlibspec( utility::pointer::dynamic_pointer_cast< rotamers::PDBRotamerLibrarySpecification const >( res_type.rotamer_library_specification() ) );
	if ( pdb_rotlibspec ) {
		libname = pdb_rotlibspec->pdb_rotamers_file();
	}
	stmt.bind(11,libname);
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
		"\tversion,\n"
		"\tname3,\n"
		"\tname1,\n"
		"\taa,\n"
		"\tlower_connect,\n"
		"\tupper_connect,\n"
		"\tnbr_atom,\n"
		"\tnbr_radius,\n"
		"\trotamer_library\n"
		"FROM residue_type\n"
		"WHERE residue_type_set_name = ? AND name = ?;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(residue_type_statement,db_session));
	stmt.bind(1,residue_type_set_name);
	stmt.bind(2,residue_type_name);
	cppdb::result res(basic::database::safely_read_from_database(stmt));

	if ( !res.next() ) {
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
	if ( version != version_ ) {
		utility_exit_with_message("Version mismatch between Residue Database and Executable");
	}

	res_type.name(residue_type_name);
	res_type.name3(name3);
	res_type.name1(name1[0]);
	res_type.aa(name_from_aa(static_cast<AA>(aa)));

	if ( lower_connect > 0 ) {
		res_type.set_lower_connect_atom(get_atom_name_from_database_atom_index(residue_type_name,lower_connect,db_session));
	}

	if ( upper_connect > 0 ) {
		res_type.set_upper_connect_atom(get_atom_name_from_database_atom_index(residue_type_name,lower_connect,db_session));
	}

	res_type.nbr_atom(get_atom_name_from_database_atom_index(residue_type_name,nbr_atom,db_session));
	res_type.nbr_radius(nbr_radius);

	if ( rotamer_library != "" ) {
		rotamers::PDBRotamerLibrarySpecificationCOP pdb_rotlibspec( utility::pointer::dynamic_pointer_cast< rotamers::PDBRotamerLibrarySpecification const >( res_type.rotamer_library_specification() ) );
		// If we have an existing library which isn't a PDB rotamers one
		if ( res_type.rotamer_library_specification() && ! pdb_rotlibspec ) {
			utility_exit_with_message("Cannot set PDB rotamer library name, restype already has a " + res_type.rotamer_library_specification()->keyname() + " library." );
		}
		res_type.rotamer_library_specification( rotamers::PDBRotamerLibrarySpecificationOP( new rotamers::PDBRotamerLibrarySpecification( rotamer_library ) ) );
	}

}

utility::vector1<std::string> ResidueDatabaseIO::get_all_residues_in_database(utility::sql_database::sessionOP db_session)
{

	utility::vector1<std::string> residue_names;

	std::string residue_name_statement =
		"SELECT\n"
		"\tname\n"
		"FROM residue_type;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(residue_name_statement,db_session));
	cppdb::result res(basic::database::safely_read_from_database(stmt));
	while ( res.next() )
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

	std::string statement_string = "INSERT INTO residue_type_atom (residue_type_set_name, residue_type_name, atom_index, atom_name, atom_type_name, mm_atom_type_name, charge, is_backbone) VALUES (?,?,?,?,?,?,?,?);";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	// AtomTypeSet?
	for ( Size i=1; i <= res_type.natoms(); ++i ) {

		stmt.bind(1,residue_type_set_name);
		stmt.bind(2,res_type.name());
		stmt.bind(3,i);
		stmt.bind(4,res_type.atom_name(i));
		stmt.bind(5,res_type.atom_type(i).atom_type_name());
		stmt.bind(6,res_type.atom(i).mm_name());
		stmt.bind(7,res_type.atom(i).charge());
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

	while ( res.next() )
			{
		core::Size atom_index;
		std::string atom_name, atom_type_name,mm_atom_type_name;
		core::Real charge;

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

	std::string statement_string = "INSERT INTO residue_type_bond (residue_type_set_name, residue_type_name, atom1, atom2, bond_type) VALUES (?,?,?,?,?);";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	for ( Size atm=1; atm <= res_type.natoms(); ++atm ) {
		AtomIndices const & neighbors(res_type.bonded_neighbor(atm));
		utility::vector1<BondName> const & types(res_type.bonded_neighbor_types(atm));
		for ( Size nbr=1; nbr <= neighbors.size(); ++nbr ) {
			if ( atm >= neighbors[nbr] ) continue;

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
		"\tatom1,\n"
		"\tatom2,\n"
		"\tbond_type\n"
		"FROM residue_type_bond\n"
		"WHERE residue_type_set_name = ? AND residue_type_name = ?;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(bond_statement_string,db_session));
	stmt.bind(1,residue_type_set_name);
	stmt.bind(2,residue_type_name);
	cppdb::result res(basic::database::safely_read_from_database(stmt));
	while ( res.next() )
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

	std::string statement_string = "INSERT INTO residue_type_cut_bond (residue_type_set_name, residue_type_name, atom1, atom2) VALUES (?,?,?,?);";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));

	for ( Size i=1; i <= res_type.natoms(); ++i ) {
		BOOST_FOREACH ( core::Size const j, res_type.cut_bond_neighbor(i) ) {
			if ( i>=j ) continue;


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
		"\tatom1,\n"
		"\tatom2\n"
		"FROM residue_type_cut_bond\n"
		"WHERE residue_type_set_name = ? AND residue_type_name = ?;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(cut_bond_statement_string,db_session));
	stmt.bind(1,residue_type_set_name);
	stmt.bind(2,residue_type_name);
	cppdb::result res(basic::database::safely_read_from_database(stmt));
	while ( res.next() )
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

	std::string statement_string = "INSERT INTO residue_type_chi (residue_type_set_name, residue_type_name, atom1, atom2, atom3, atom4, chino) VALUES (?,?,?,?,?,?,?);";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));

	for ( Size i=1; i <= res_type.nchi(); ++i ) {
		AtomIndices const & chi_atoms(res_type.chi_atoms(i));

		stmt.bind(1,residue_type_set_name);
		stmt.bind(2,res_type.name());
		stmt.bind(3,chi_atoms[1]);
		stmt.bind(4,chi_atoms[2]);
		stmt.bind(5,chi_atoms[3]);
		stmt.bind(6,chi_atoms[4]);
		stmt.bind(7,i);
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
		"\tchino,\n"
		"\tatom1,\n"
		"\tatom2,\n"
		"\tatom3,\n"
		"\tatom4\n"
		"FROM residue_type_chi\n"
		"WHERE residue_type_set_name = ? AND residue_type_name = ?;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(chi_statement_string,db_session));
	stmt.bind(1,residue_type_set_name);
	stmt.bind(2,residue_type_name);
	cppdb::result res(basic::database::safely_read_from_database(stmt));
	while ( res.next() )
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

	std::string statement_string = "INSERT INTO residue_type_chi_rotamer (residue_type_set_name, residue_type_name, chino, mean, sdev) VALUES (?,?,?,?,?);";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	for ( Size chi=1; chi <= res_type.nchi(); ++chi ) {
		std::pair<Real, Real> mean_sdev;
		BOOST_FOREACH ( mean_sdev, res_type.chi_rotamers(chi) ) {

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
		"\tchino,\n"
		"\tmean,\n"
		"\tsdev\n"
		"FROM residue_type_chi_rotamer\n"
		"WHERE residue_type_set_name = ? AND residue_type_name = ?;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(chi_rotamer_statement_string,db_session));
	stmt.bind(1,residue_type_set_name);
	stmt.bind(2,residue_type_name);
	cppdb::result res(basic::database::safely_read_from_database(stmt));
	while ( res.next() )
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

	std::string statement_string = "INSERT INTO residue_type_proton_chi (residue_type_set_name, residue_type_name, chino, sample, is_extra) VALUES (?,?,?,?,?);";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	for ( Size proton_chi=1; proton_chi <= res_type.n_proton_chi(); ++proton_chi ) {
		Size const chi(res_type.proton_chi_2_chi(proton_chi));
		BOOST_FOREACH ( Real const sample, res_type.proton_chi_samples(proton_chi) ) {

			stmt.bind(1,residue_type_set_name);
			stmt.bind(2,res_type.name());
			stmt.bind(3,chi);
			stmt.bind(4,sample);
			stmt.bind(5,false);
			basic::database::safely_write_to_database(stmt);

			BOOST_FOREACH ( Real const extra_sample, res_type.proton_chi_extra_samples(proton_chi) ) {

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
	core::chemical::ResidueType &,
	utility::sql_database::sessionOP db_session)
{
	std::string proton_chi_statement_string =
		"SELECT\n"
		"\tchino,\n"
		"\tsample,\n"
		"\tis_extra\n"
		"FROM residue_type_proton_chi\n"
		"WHERE residue_type_set_name = ? AND residue_type_name = ?;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(proton_chi_statement_string,db_session));
	stmt.bind(1,residue_type_set_name);
	stmt.bind(2,residue_type_name);
	cppdb::result res(basic::database::safely_read_from_database(stmt));

	std::map<core::Size,utility::vector1< core::Real > > samples;
	std::map<core::Size,utility::vector1< core::Real > > extra_samples;

	while ( res.next() )
			{
		core::Size chino, is_extra;
		core::Real sample;

		res >> chino >> sample >> is_extra;
		if ( is_extra ) {
			if ( extra_samples.find(chino) == extra_samples.end() ) {
				utility::vector1< core::Real> extra_sample_vect;
				extra_sample_vect.push_back(sample);
				extra_samples.insert(std::make_pair(chino,extra_sample_vect));
			} else {
				extra_samples[chino].push_back(sample);
			}
		} else {
			//we've never seen this chino, add it to the map
			if ( samples.find(chino) == samples.end() ) {
				utility::vector1< core::Real> sample_vect;
				sample_vect.push_back(sample);
				samples.insert(std::make_pair(chino,sample_vect));
			} else {
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
	BOOST_FOREACH ( std::string const & property, res_type.properties().get_list_of_properties() ) {

		cppdb::statement stmt = (*db_session)
			<< "INSERT INTO residue_type_property (residue_type_set_name, residue_type_name, property) VALUES (?,?,?);"
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
		"\tproperty\n"
		"FROM residue_type_property\n"
		"WHERE residue_type_set_name = ? AND residue_type_name = ?;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(residue_property_statement_string,db_session));
	stmt.bind(1,residue_type_set_name);
	stmt.bind(2,residue_type_name);
	cppdb::result res(basic::database::safely_read_from_database(stmt));
	while ( res.next() )
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
	utility::sql_database::sessionOP db_session )
{
	BOOST_FOREACH ( std::string const & variant_type, res_type.properties().get_list_of_variants() ) {
		cppdb::statement stmt = (*db_session)
			<< "INSERT INTO residue_type_variant_type (residue_type_set_name, residue_type_name, variant_type) VALUES (?,?,?);"
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
		"\tvariant_type\n"
		"FROM residue_type_variant_type\n"
		"WHERE residue_type_set_name = ? AND residue_type_name = ?;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(residue_variant_statement_string,db_session));
	stmt.bind(1,residue_type_set_name);
	stmt.bind(2,residue_type_name);
	cppdb::result res(basic::database::safely_read_from_database(stmt));
	while ( res.next() )
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
	std::string statement_string = "INSERT INTO residue_type_icoor (residue_type_set_name, residue_type_name, child_atom, icoor_sequence, phi, theta, distance, parent_atom, angle_atom, torsion_atom) VALUES (?,?,?,?,?,?,?,?,?,?);";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));

	for ( Size i=1; i <= res_type.natoms(); ++i ) {

		AtomICoor atom_icoor(res_type.icoor(i));
		stmt.bind(1,residue_type_set_name);
		stmt.bind(2,res_type.name());
		stmt.bind(3,res_type.atom_name(i));
		// This was atom_icoor.index(), but that no longer exists.
		// Retain position for backward compatability.
		stmt.bind(4,i);
		stmt.bind(5,atom_icoor.phi());
		stmt.bind(6,atom_icoor.theta());
		stmt.bind(7,atom_icoor.d());
		//We need to do this because the AtomICoor doesn't store atom names, just numbers, but it takes atom names as its constructor
		if ( atom_icoor.stub_atom1().type() == ICoorAtomID::INTERNAL ) {
			stmt.bind(8,res_type.atom_name(atom_icoor.stub_atom1().atomno()));
		} else if ( atom_icoor.stub_atom1().type() == ICoorAtomID::CONNECT ) {
			stmt.bind(8,"CONNECT");
		} else if ( atom_icoor.stub_atom2().type() == ICoorAtomID::POLYMER_LOWER ) {
			stmt.bind(8,"LOWER");
		} else {
			//should be POLYMER_UPPER
			stmt.bind(8,"UPPER");
		}

		if ( atom_icoor.stub_atom2().type() == ICoorAtomID::INTERNAL ) {
			stmt.bind(9,res_type.atom_name(atom_icoor.stub_atom2().atomno()));
		} else if ( atom_icoor.stub_atom1().type() == ICoorAtomID::CONNECT ) {
			stmt.bind(9,"CONNECT");
		} else if ( atom_icoor.stub_atom2().type() == ICoorAtomID::POLYMER_LOWER ) {
			stmt.bind(9,"LOWER");
		} else {
			//should be POLYMER_UPPER
			stmt.bind(9,"UPPER");
		}

		if ( atom_icoor.stub_atom3().type() == ICoorAtomID::INTERNAL ) {
			stmt.bind(10,res_type.atom_name(atom_icoor.stub_atom3().atomno()));
		} else if ( atom_icoor.stub_atom1().type() == ICoorAtomID::CONNECT ) {
			stmt.bind(10,"CONNECT");
		} else if ( atom_icoor.stub_atom2().type() == ICoorAtomID::POLYMER_LOWER ) {
			stmt.bind(10,"LOWER");
		} else {
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
		"\tchild_atom,\n"
		"\tphi,\n"
		"\ttheta,\n"
		"\tdistance,\n"
		"\tparent_atom,\n"
		"\tangle_atom,\n"
		"\ttorsion_atom\n"
		"FROM residue_type_icoor\n"
		"WHERE residue_type_set_name = ? AND residue_type_name = ?\n"
		"ORDER BY icoor_sequence;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(residue_variant_statement_string,db_session));
	stmt.bind(1,residue_type_set_name);
	stmt.bind(2,residue_type_name);

	cppdb::result res(basic::database::safely_read_from_database(stmt));
	while ( res.next() )
			{
		std::string child_atom_name, parent_atom_name,angle_atom_name,torsion_atom_name;
		core::Real phi, theta, distance;
		res >> child_atom_name >> phi >> theta >> distance >> parent_atom_name >>angle_atom_name >>torsion_atom_name;

		// set icoor
		// (Atom base gets set in set_icoor)
		res_type.set_icoor( child_atom_name, phi, theta, distance, parent_atom_name, angle_atom_name, torsion_atom_name );
	}

	// We do this at the end of the loop to make sure all the dependant icoors are present.
	core::chemical::fill_ideal_xyz_from_icoor( res_type, res_type.graph() );

}


}
}
