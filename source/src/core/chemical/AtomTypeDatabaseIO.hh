// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/chemical/AtomTypeDatabaseIO.hh
/// @author Matthew O'Meara

#ifndef INCLUDED_core_chemical_AtomTypeDatabaseIO_hh
#define INCLUDED_core_chemical_AtomTypeDatabaseIO_hh

#include <utility/pointer/ReferenceCount.hh>
//#include <core/chemical/ResidueType.fwd.hh>
//#include <core/chemical/ResidueDatabaseIO.fwd.hh>

// needed for header only build
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>
#include <utility/vector1.hh>

#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>

//#include <map>

namespace core {
namespace chemical {


class AtomTypeDatabaseIO : public utility::pointer::ReferenceCount {
public:

	AtomTypeDatabaseIO();

	~AtomTypeDatabaseIO() override;

public: // store in a database information associated with an atom type set

	/// @brief generate the table schemas and write them to the database
	void
	write_schema_to_db(
		utility::sql_database::sessionOP db_session) const;

private:
	void
	write_atom_types_table_schema(
		utility::sql_database::sessionOP db_session) const;

	void
	write_atom_type_property_values_table_schema(
		utility::sql_database::sessionOP db_session) const;

	void
	write_atom_type_properties_table_schema(
		utility::sql_database::sessionOP db_session) const;

	void
	write_atom_type_extra_parameters_table_schema(
		utility::sql_database::sessionOP db_session) const;


public:
	/// @brief write the schema
	void
	initialize(
		utility::sql_database::sessionOP db_session) const;

	void
	write_atom_type_set_to_database(
		chemical::AtomTypeSet const & atom_type_set,
		utility::sql_database::sessionOP db_session) const;

	utility::vector1<std::string>
	get_all_atom_types_in_database(
		utility::sql_database::sessionOP db_session) const;

private:

	void
	write_atom_type_table(
		std::string const & atom_type_set_name,
		AtomType const & atom_type,
		utility::sql_database::sessionOP db_session) const;

	void
	write_atom_type_properties_table(
		std::string const & atom_type_set_name,
		AtomType const & atom_type,
		utility::sql_database::sessionOP db_session) const;

	void
	write_atom_type_extra_parameters_table(
		AtomTypeSet const & atom_type_set,
		AtomType const & atom_type,
		utility::sql_database::sessionOP db_session) const;

};


} // namespace
} // namespace

#endif // include guard
