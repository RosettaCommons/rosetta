// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ResidueTypesFeatures.cc
/// @brief  report ResidueTypes to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/ResidueTypesFeatures.hh>

// Project Headers
#include <basic/Tracer.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <basic/database/sql_utils.hh>

// Utility Headers
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>

// External Headers
#include <cppdb/frontend.h>

// Boost Headers
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// C++ Headers
#include <set>
#include <sstream>

namespace protocols{
namespace features{

using std::string;
using std::stringstream;
using std::endl;
using std::pair;
using std::set;
using core::Size;
using core::Real;
using core::pose::Pose;
using core::chemical::BondName;
using core::chemical::AtomIndices;
using core::chemical::ResidueType;
using utility::sql_database::sessionOP;
using utility::vector1;
using cppdb::statement;
using cppdb::result;
using basic::Tracer;

static Tracer TR("protocols.features.ResidueTypesFeatures");

ResidueTypesFeatures::ResidueTypesFeatures() :
	version_(0.1) // should match version string in residue type parameter sets
{}

ResidueTypesFeatures::~ResidueTypesFeatures() {}

string
ResidueTypesFeatures::type_name() const { return "ResidueTypesFeatures"; }

string
ResidueTypesFeatures::schema() const {
	// NOTE: To support building feature databases in parallel, the
	// ResidueTypeSet and ResidueType objects must be identified by
	// their names rather then assigning them a unique id.
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
		"	PRIMARY KEY(residue_type_set_name, residue_type_name, atom1, atom2));\n"
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
		"	PRIMARY KEY(residue_type_set_name, residue_type_name, variant_type));\n";
}

Size
ResidueTypesFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size const struct_id,
	sessionOP db_session
){

	// Get a set of the unique residue types that are used in this pose
	set< ResidueType const * > res_types;
	for(Size i=1; i <= pose.n_residue(); ++i){
		if(!relevant_residues[i]) continue;
		res_types.insert(&pose.residue_type(i));
	}

	foreach(ResidueType const * res_type, res_types){
		string const & residue_type_set_name(res_type->residue_type_set().name());

		// Is this residue type already in the database?
		result res;
		while(true)
		{
			try
			{
				res = (*db_session)
					<< "SELECT * FROM residue_type\n"
					"WHERE residue_type_set_name = ? AND name = ?;"
					<< residue_type_set_name << res_type->name();
				break;
			}catch(cppdb::cppdb_error &)
			{
				#ifndef WIN_PYROSETTA
					usleep(10);
				#endif
				continue;
			}
		}
		if(res.next()) continue;


		report_residue_type(residue_type_set_name, *res_type, db_session);
		report_residue_type_atom(residue_type_set_name, *res_type, db_session);
		report_residue_type_bond(residue_type_set_name, *res_type, db_session);
		report_residue_type_cut_bond(residue_type_set_name, *res_type, db_session);
		report_residue_type_chi(residue_type_set_name, *res_type, db_session);
		report_residue_type_chi_rotamer(residue_type_set_name, *res_type, db_session);
		report_residue_type_proton_chi(residue_type_set_name, *res_type, db_session);
		report_residue_type_properties(residue_type_set_name, *res_type, db_session);
		report_residue_type_variant(residue_type_set_name, *res_type, db_session);
	}
	return 0;
}

void
ResidueTypesFeatures::report_residue_type(
	string const & residue_type_set_name,
	ResidueType const & res_type,
	sessionOP db_session
) const {

	stringstream name1; name1 << res_type.name1();

	int lower_terminus(-1), upper_terminus(-1);
	if(res_type.is_polymer()){
		if(!res_type.is_lower_terminus()) lower_terminus = res_type.lower_connect_atom();
		if(!res_type.is_upper_terminus()) upper_terminus = res_type.upper_connect_atom();
	}


	statement stmt = (*db_session)
		<< "INSERT INTO residue_type VALUES (?,?,?,?,?,?,?,?,?,?);"
		<< residue_type_set_name
		<< version_
		<< res_type.name()
		<< res_type.name3()
		<< name1.str()
		<< res_type.aa()
		<< lower_terminus
		<< upper_terminus
		<< res_type.nbr_atom()
		<< res_type.nbr_radius();
	basic::database::safely_write_to_database(stmt);

}

void
ResidueTypesFeatures::report_residue_type_atom(
	string const & residue_type_set_name,
	ResidueType const & res_type,
	sessionOP db_session
) const {

	// AtomTypeSet?
	for(Size i=1; i <= res_type.natoms(); ++i){

		statement stmt = (*db_session)
			<< "INSERT INTO residue_type_atom VALUES (?,?,?,?,?,?,?,?);"
			<< residue_type_set_name
			<< res_type.name()
			<< i
			<< res_type.atom_name(i)
			<< res_type.atom_type(i).atom_type_name()
			<< res_type.mm_atom_name(i)
			<< res_type.atomic_charge(i)
			<< res_type.atom_is_backbone(i);
		basic::database::safely_write_to_database(stmt);

	}
}

void
ResidueTypesFeatures::report_residue_type_bond(
	string const & residue_type_set_name,
	ResidueType const & res_type,
	sessionOP db_session
) const {

	for(Size atm=1; atm <= res_type.natoms(); ++atm){
		AtomIndices const & neighbors(res_type.bonded_neighbor(atm));
		vector1<BondName> const & types(res_type.bonded_neighbor_types(atm));
		for(Size nbr=1; nbr <= neighbors.size(); ++nbr){
			if(atm >= neighbors[nbr]) continue;

			statement stmt = (*db_session)
				<< "INSERT INTO residue_type_bond VALUES (?,?,?,?,?);"
				<< residue_type_set_name
				<< res_type.name()
				<< atm
				<< neighbors[nbr]
				<< types[nbr];
			basic::database::safely_write_to_database(stmt);
		}
	}
}

void
ResidueTypesFeatures::report_residue_type_cut_bond(
	string const & residue_type_set_name,
	ResidueType const & res_type,
	sessionOP db_session
) const {

	for(Size i=1; i <= res_type.natoms(); ++i){
		foreach(Size const j, res_type.cut_bond_neighbor(i)){
			if(i>=j) continue;
			statement stmt = (*db_session)
				<< "INSERT INTO residue_type_cut_bond VALUES (?,?,?,?);"
				<< residue_type_set_name
				<< res_type.name()
				<< i
				<< j;
			basic::database::safely_write_to_database(stmt);
		}
	}
}


void
ResidueTypesFeatures::report_residue_type_chi(
	string const & residue_type_set_name,
	ResidueType const & res_type,
	sessionOP db_session
) const {
	for(Size i=1; i <= res_type.nchi(); ++i){
		AtomIndices const & chi_atoms(res_type.chi_atoms(i));

		statement stmt = (*db_session)
			<< "INSERT INTO residue_type_chi VALUES (?,?,?,?,?,?,?);"
			<< residue_type_set_name
			<< res_type.name()
			<< i
			<< chi_atoms[1]
			<< chi_atoms[2]
			<< chi_atoms[3]
			<< chi_atoms[4];
		basic::database::safely_write_to_database(stmt);
	}

}

void
ResidueTypesFeatures::report_residue_type_chi_rotamer(
	string const & residue_type_set_name,
	ResidueType const & res_type,
	sessionOP db_session
) const {
	for(Size chi=1; chi <= res_type.nchi(); ++chi){
		pair<Real, Real> mean_sdev;
		foreach(mean_sdev, res_type.chi_rotamers(chi)){

			statement stmt = (*db_session)
				<< "INSERT INTO residue_type_chi_rotamer VALUES (?,?,?,?,?);"
				<< residue_type_set_name
				<< res_type.name()
				<< chi
				<< mean_sdev.first
				<< mean_sdev.second;
			basic::database::safely_write_to_database(stmt);
		}
	}
}

void
ResidueTypesFeatures::report_residue_type_proton_chi(
	string const & residue_type_set_name,
	ResidueType const & res_type,
	sessionOP db_session
) const {
	for(Size proton_chi=1; proton_chi <= res_type.n_proton_chi(); ++proton_chi){
		Size const chi(res_type.proton_chi_2_chi(proton_chi));
		foreach(Real const sample, res_type.proton_chi_samples(proton_chi)){

			statement stmt_sample = (*db_session)
				<< "INSERT INTO residue_type_proton_chi VALUES (?,?,?,?,?);"
				<< residue_type_set_name
				<< res_type.name()
				<< chi
				<< sample
				<< false;
			basic::database::safely_write_to_database(stmt_sample);


			foreach(Real const extra_sample, res_type.proton_chi_extra_samples(proton_chi)){

				statement stmt_neg_extra = (*db_session)
					<< "INSERT INTO residue_type_proton_chi VALUES (?,?,?,?,?);"
					<< residue_type_set_name
					<< res_type.name()
					<< chi
					<< sample - extra_sample
					<< true;
				basic::database::safely_write_to_database(stmt_neg_extra);

				statement stmt_pos_extra = (*db_session)
					<< "INSERT INTO residue_type_proton_chi VALUES (?,?,?,?,?);"
					<< residue_type_set_name
					<< res_type.name()
					<< chi
					<< sample + extra_sample
					<< true;
				basic::database::safely_write_to_database(stmt_pos_extra);

			}
		}
	}
}

void
ResidueTypesFeatures::report_residue_type_properties(
	string const & residue_type_set_name,
	ResidueType const & res_type,
	sessionOP db_session
) const {
	foreach(string const & property, res_type.properties()){

				statement stmt = (*db_session)
					<< "INSERT INTO residue_type_property VALUES (?,?,?);"
					<< residue_type_set_name
					<< res_type.name()
					<< property;
				basic::database::safely_write_to_database(stmt);

	}
}

void
ResidueTypesFeatures::report_residue_type_variant(
	string const & residue_type_set_name,
	ResidueType const & res_type,
	sessionOP db_session
) const {
	foreach(string const & variant_type, res_type.variant_types()){

		statement stmt = (*db_session)
			<< "INSERT INTO residue_type_variant_type VALUES (?,?,?);"
			<< residue_type_set_name
			<< res_type.name()
			<< variant_type;
		basic::database::safely_write_to_database(stmt);

	}
}

} // namesapce
} // namespace
