// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ResidueConformationFeatures.cc
/// @brief  report idealized torsional DOFs Statistics Scientific Benchmark
/// @author Matthew O'Meara
/// @author Sam DeLuca

// Unit Headers
#include <protocols/features/ResidueConformationFeatures.hh>

//External
#include <boost/uuid/uuid.hpp>

// Project Headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/id/AtomID.hh>

// Utility Headers
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

//Basic Headers
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/Constraint.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

// External Headers
#include <cppdb/frontend.h>
#include <boost/uuid/uuid_io.hpp>

// C++ Headers
#include <cmath>
#include <sstream>

namespace protocols{
namespace features{

using std::string;
using core::Size;
using core::Real;
using core::pose::Pose;
using core::conformation::Residue;
using core::chemical::num_canonical_aas;
using basic::database::table_exists;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;

string
ResidueConformationFeatures::type_name() const {
	return "ResidueConformationFeatures";
}

void
ResidueConformationFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const{
	using namespace basic::database::schema_generator;

	//******nonprotein_residue_conformation******//
	Column struct_id("struct_id",DbUUID(), false);
	Column seqpos("seqpos",DbInteger(), false);
	Column phi("phi",DbDouble(), false);
	Column psi("psi",DbDouble(), false);
	Column omega("omega",DbDouble(), false);

	utility::vector1<Column> non_prot_res_pkeys;
	non_prot_res_pkeys.push_back(struct_id);
	non_prot_res_pkeys.push_back(seqpos);



	utility::vector1<Column> fkey_cols;
	fkey_cols.push_back(struct_id);
	fkey_cols.push_back(seqpos);

	utility::vector1<std::string> fkey_reference_cols;
	fkey_reference_cols.push_back("struct_id");
	fkey_reference_cols.push_back("resNum");

	Schema nonprotein_residue_conformation("nonprotein_residue_conformation", PrimaryKey(non_prot_res_pkeys));
	nonprotein_residue_conformation.add_column(struct_id);
	nonprotein_residue_conformation.add_column(seqpos);
	nonprotein_residue_conformation.add_column(phi);
	nonprotein_residue_conformation.add_column(psi);
	nonprotein_residue_conformation.add_column(omega);
	nonprotein_residue_conformation.add_foreign_key(ForeignKey(fkey_cols, "residues", fkey_reference_cols, true));

	nonprotein_residue_conformation.write(db_session);

	//******nonprotein_residue_angles******//
	Column chinum("chinum",DbInteger(), false);
	Column chiangle("chiangle",DbDouble(), false);

	utility::vector1<Column> non_prot_res_angle_keys;
	non_prot_res_angle_keys.push_back(struct_id);
	non_prot_res_angle_keys.push_back(seqpos);
	non_prot_res_angle_keys.push_back(chinum);

	Schema nonprotein_residue_angles("nonprotein_residue_angles", PrimaryKey(non_prot_res_angle_keys));
	nonprotein_residue_angles.add_column(struct_id);
	nonprotein_residue_angles.add_column(seqpos);
	nonprotein_residue_angles.add_column(chinum);
	nonprotein_residue_angles.add_column(chiangle);
	nonprotein_residue_angles.add_foreign_key(ForeignKey(fkey_cols, "residues", fkey_reference_cols, true));

	nonprotein_residue_angles.write(db_session);

	//******residue_atom_coords******//
	Column atomno("atomno",DbInteger(), false);
	Column x("x",DbDouble(), false);
	Column y("y",DbDouble(), false);
	Column z("z",DbDouble(), false);

	utility::vector1<Column> res_atm_coords_pkeys;
	res_atm_coords_pkeys.push_back(struct_id);
	res_atm_coords_pkeys.push_back(seqpos);
	res_atm_coords_pkeys.push_back(atomno);

	Schema residue_atom_coords("residue_atom_coords", PrimaryKey(res_atm_coords_pkeys));
	residue_atom_coords.add_column(struct_id);
	residue_atom_coords.add_column(seqpos);
	residue_atom_coords.add_column(atomno);
	residue_atom_coords.add_column(x);
	residue_atom_coords.add_column(y);
	residue_atom_coords.add_column(z);
	residue_atom_coords.add_foreign_key(ForeignKey(fkey_cols, "residues", fkey_reference_cols, true));

	residue_atom_coords.write(db_session);

//	if(db_mode == "sqlite3")
//	{
//		return
//			"CREATE TABLE IF NOT EXISTS nonprotein_residue_conformation (\n"
//			"	struct_id BLOB,\n"
//			"	seqpos INTEGER,\n"
//			"	phi REAL,\n"
//			"	psi REAL,\n"
//			"	omega REAL,\n"
//			"	FOREIGN KEY (struct_id, seqpos)\n"
//			"		REFERENCES residues (struct_id, resNum)\n"
//			"		DEFERRABLE INITIALLY DEFERRED,\n"
//			"	PRIMARY KEY (struct_id, seqpos));\n"
//			"\n"
//			"CREATE TABLE IF NOT EXISTS nonprotein_residue_angles (\n"
//			"	struct_id BLOB,\n"
//			"	seqpos INTEGER,\n"
//			"	chinum INTEGER,\n"
//			"	chiangle REAL,\n"
//			"	FOREIGN KEY (struct_id, seqpos)\n"
//			"		REFERENCES residues (struct_id, resNum)\n"
//			"		DEFERRABLE INITIALLY DEFERRED,\n"
//			"	PRIMARY KEY (struct_id, seqpos, chinum));\n"
//			"\n"
//			"CREATE TABLE IF NOT EXISTS residue_atom_coords (\n"
//			"	struct_id BLOB,\n"
//			"	seqpos INTEGER,\n"
//			"	atomno INTEGER,\n"
//			"	x REAL,\n"
//			"	y REAL,\n"
//			"	z REAL,\n"
//			"	FOREIGN KEY (struct_id, seqpos)\n"
//			"		REFERENCES residues (struct_id, resNum)\n"
//			"		DEFERRABLE INITIALLY DEFERRED,\n"
//			"	PRIMARY KEY (struct_id, seqpos, atomno));\n";
//	}else if(db_mode == "mysql")
//	{
//		return
//			"CREATE TABLE IF NOT EXISTS nonprotein_residue_conformation (\n"
//			"	struct_id BINARY(16),\n"
//			"	seqpos INTEGER,\n"
//			"	phi DOUBLE,\n"
//			"	psi DOUBLE,\n"
//			"	omega DOUBLE,\n"
//			"	FOREIGN KEY (struct_id, seqpos) REFERENCES residues (struct_id, resNum),"
//			"	PRIMARY KEY (struct_id, seqpos));\n"
//			"\n"
//			"CREATE TABLE IF NOT EXISTS nonprotein_residue_angles (\n"
//			"	struct_id BINARY(16),\n"
//			"	seqpos INTEGER,\n"
//			"	chinum INTEGER,\n"
//			"	chiangle DOUBLE,\n"
//			"	FOREIGN KEY (struct_id, seqpos) REFERENCES residues (struct_id, resNum),"
//			"	PRIMARY KEY (struct_id, seqpos, chinum));\n"
//			"\n"
//			"CREATE TABLE IF NOT EXISTS residue_atom_coords (\n"
//			"	struct_id BINARY(16),\n"
//			"	seqpos INTEGER,\n"
//			"	atomno INTEGER,\n"
//			"	x DOUBLE,\n"
//			"	y DOUBLE,\n"
//			"	z DOUBLE,\n"
//			"	FOREIGN KEY (struct_id, seqpos) REFERENCES residues (struct_id, resNum),"
//			"	PRIMARY KEY (struct_id, seqpos, atomno));";
//	}else
//	{
//		return "";
//	}

}

utility::vector1<std::string>
ResidueConformationFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	return dependencies;
}

Size
ResidueConformationFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	boost::uuids::uuid struct_id,
	sessionOP db_session
){
	//check to see if this structure is ideal

	bool ideal = true;
	if(!basic::options::option[basic::options::OptionKeys::out::file::force_nonideal_structure]())
	{
		core::conformation::Conformation const & conformation(pose.conformation());
		for(core::Size resn=1; resn <= pose.n_residue();++resn){
			if(!relevant_residues[resn]) continue;
			bool residue_status(core::conformation::is_ideal_position(resn,conformation));
			if(!residue_status){
				ideal = false;
				break;
			}
		}
	}else
	{
		ideal = false;
	}

	std::string conformation_string = "INSERT INTO nonprotein_residue_conformation (struct_id, seqpos, phi, psi, omega) VALUES (?,?,?,?,?)";
	std::string angle_string = "INSERT INTO nonprotein_residue_angles (struct_id, seqpos, chinum, chiangle) VALUES (?,?,?,?)";
	std::string coords_string = "INSERT INTO residue_atom_coords (struct_id, seqpos, atomno, x, y, z) VALUES (?,?,?,?,?,?)";

	statement conformation_stmt(basic::database::safely_prepare_statement(conformation_string,db_session));
	statement angle_stmt(basic::database::safely_prepare_statement(angle_string,db_session));
	statement coords_stmt(basic::database::safely_prepare_statement(coords_string,db_session));


	//cppdb::transaction transact_guard(*db_session);
	for (Size i = 1; i <= pose.total_residue(); ++i) {
		if(!relevant_residues[i]) continue;

		Residue const & resi = pose.residue(i);
		if(resi.aa() <= num_canonical_aas){
			continue;
		}
		//runtime_assert(resi.aa() <= num_canonical_aas);
		Real phi  (0.0);
		Real psi  (0.0);
		Real omega(0.0);
		if(!resi.is_ligand()){
			phi  = resi.mainchain_torsion(1);
			psi  = resi.mainchain_torsion(2);
			omega= resi.mainchain_torsion(3);
		}

		conformation_stmt.bind(1,struct_id);
		conformation_stmt.bind(2,i);
		conformation_stmt.bind(3,phi);
		conformation_stmt.bind(4,psi);
		conformation_stmt.bind(5,omega);
		basic::database::safely_write_to_database(conformation_stmt);

		for(core::Size chi_num = 1; chi_num <= resi.nchi();++chi_num){
			core::Real chi_angle = resi.chi(chi_num);

			angle_stmt.bind(1,struct_id);
			angle_stmt.bind(2,i);
			angle_stmt.bind(3,chi_num);
			angle_stmt.bind(4,chi_angle);
			basic::database::safely_write_to_database(angle_stmt);

		}
		if(!ideal || resi.is_ligand()){ // always store coords for a ligand
			for(Size atom = 1; atom <= resi.natoms(); ++atom){
				core::Vector coords = resi.xyz(atom);

				coords_stmt.bind(1,struct_id);
				coords_stmt.bind(2,i);
				coords_stmt.bind(3,atom);
				coords_stmt.bind(4,coords.x());
				coords_stmt.bind(5,coords.y());
				coords_stmt.bind(6,coords.z());
				basic::database::safely_write_to_database(coords_stmt);

			}
		}
	}
	//transact_guard.commit();
	return 0;
}

void
ResidueConformationFeatures::delete_record(
	boost::uuids::uuid struct_id,
	sessionOP db_session
){

	statement conformation_stmt(basic::database::safely_prepare_statement("DELETE FROM nonprotein_residue_conformation WHERE struct_id = ?;\n",db_session));
	conformation_stmt.bind(1,struct_id);
	basic::database::safely_write_to_database(conformation_stmt);
	statement angle_stmt(basic::database::safely_prepare_statement("DELETE FROM nonprotein_residue_angles WHERE struct_id = ?;\n",db_session));
	angle_stmt.bind(1,struct_id);
	basic::database::safely_write_to_database(angle_stmt);
	statement coords_stmt(basic::database::safely_prepare_statement("DELETE FROM residue_atom_coords WHERE struct_id = ?;",db_session));
	coords_stmt.bind(1,struct_id);
	basic::database::safely_write_to_database(coords_stmt);

}

void
ResidueConformationFeatures::load_into_pose(
	sessionOP db_session,
	boost::uuids::uuid struct_id,
	Pose & pose
){
	load_conformation(db_session, struct_id, pose);
}

void
ResidueConformationFeatures::load_conformation(
	sessionOP db_session,
	boost::uuids::uuid struct_id,
	Pose & pose
){
	if(!table_exists(db_session, "nonprotein_residue_conformation")) return;
	if(!table_exists(db_session, "nonprotein_residue_angles")) return;

	if(pose.is_fullatom()){

		std::string protein_string =
			"SELECT\n"
			"	seqpos,\n"
			"	phi,\n"
			"	psi,\n"
			"	omega\n"
			"FROM\n"
			"	nonprotein_residue_conformation\n"
			"WHERE\n"
			"	nonprotein_residue_conformation.struct_id=?;";

		statement protein_stmt(basic::database::safely_prepare_statement(protein_string,db_session));
		protein_stmt.bind(1,struct_id);

		std::string conformation_string =
			"SELECT\n"
			"	seqpos,\n"
			"	chinum,\n"
			"	chiangle\n"
			"FROM\n"
			"	nonprotein_residue_angles\n"
			"WHERE\n"
			"	nonprotein_residue_angles.struct_id=?;";

		statement conformation_stmt(basic::database::safely_prepare_statement(conformation_string,db_session));
		conformation_stmt.bind(1,struct_id);


		result protein_res(basic::database::safely_read_from_database(protein_stmt));
		while(protein_res.next()){
			Size seqpos;
			Real phi,psi,omega;
			protein_res >> seqpos >> phi >> psi >> omega;
			if (pose.residue_type(seqpos).is_protein()){
				pose.set_phi(seqpos,phi);
				pose.set_psi(seqpos,psi);
				pose.set_omega(seqpos,omega);
			}
		}
		result res_conformation(basic::database::safely_read_from_database(conformation_stmt));
		while(res_conformation.next()){
			//Size nchi(pose.residue_type(seqpos).nchi());
			Size seqpos;
			Size chinum;
			Real chiangle;
			set_coords_for_residue(db_session,struct_id,seqpos,pose);
			res_conformation >> seqpos >> chinum >> chiangle;
			pose.set_chi(chinum,seqpos,chiangle);
		}

	}else{

		std::string protein_string  =
			"SELECT\n"
			"	seqpos,\n"
			"	phi,\n"
			"	psi,\n"
			"	omega\n"
			"FROM\n"
			"	nonprotein_residue_conformation\n"
			"WHERE\n"
			"	nonprotein_residue_conformation.struct_id=?;";

		statement protein_stmt(basic::database::safely_prepare_statement(protein_string,db_session));
		protein_stmt.bind(1,struct_id);

		result protein_res(basic::database::safely_read_from_database(protein_stmt));
		while(protein_res.next()){
			Size seqpos;
			Real phi,psi,omega;
			protein_res >> seqpos >> phi >> psi >> omega;
			if (!pose.residue_type(seqpos).is_protein()){
				// WARNING why are you storing non-protein in the ProteinSilentReport?
				continue;
			}
			set_coords_for_residue(db_session,struct_id,seqpos,pose);
			pose.set_phi(seqpos,phi);
			pose.set_psi(seqpos,psi);
			pose.set_omega(seqpos,omega);
		}

	}
}

//This should be factored out into a non-member function.
void ResidueConformationFeatures::set_coords_for_residue(
		sessionOP db_session,
		boost::uuids::uuid struct_id,
		Size seqpos,
		Pose & pose
){

	std::string statement_string =
		"SELECT\n"
		"	atomno,\n"
		"	x,\n"
		"	y,\n"
		"	z\n"
		"FROM\n"
		"	residue_atom_coords\n"
		"WHERE\n"
		"residue_atom_coords.struct_id=? AND residue_atom_coords.seqpos=?;";

	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);
	stmt.bind(2,seqpos);

	result res(basic::database::safely_read_from_database(stmt));
	while(res.next()){
		Size atomno;
		Real x,y,z;
		res >> atomno >> x >> y >> z;

		core::id::AtomID atom_id(atomno,seqpos);
		core::Vector coords(x,y,z);
		pose.set_xyz(atom_id,coords);
	}


}


} // features
} // protocols
