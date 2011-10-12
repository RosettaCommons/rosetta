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

// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>

// External Headers
#include <cppdb/frontend.h>

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
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;

string
ResidueConformationFeatures::type_name() const {
	return "ProteinResidueConformationFeatures";
}

string
ResidueConformationFeatures::schema() const {
	std::string db_mode(basic::options::option[basic::options::OptionKeys::inout::database_mode]);

	if(db_mode == "sqlite3")
	{
		return
			"CREATE TABLE IF NOT EXISTS nonprotein_residue_conformation (\n"
			"	struct_id INTEGER,\n"
			"	seqpos INTEGER,\n"
			"	phi REAL,\n"
			"	psi REAL,\n"
			"	omega REAL,\n"
			"	FOREIGN KEY (struct_id)\n"
			"		REFERENCES structures (struct_id)\n"
			"		DEFERRABLE INITIALLY DEFERRED);"
			"\n"
			"CREATE TABLE IF NOT EXISTS nonprotein_residue_angles (\n"
			"	struct_id INTEGER,\n"
			"	seqpos INTEGER,\n"
			"	chinum INTEGER,\n"
			"	chiangle REAL,\n"
			"	FOREIGN KEY (struct_id)\n"
			"		REFERENCES structures (struct_id)\n"
			"		DEFERRABLE INITIALLY DEFERRED);"
			"\n"
			"CREATE TABLE IF NOT EXISTS residue_atom_coords (\n"
			"	struct_id INTEGER,\n"
			"	seqpos INTEGER,\n"
			"	atomno INTEGER,\n"
			"	x REAL,\n"
			"	y REAL,\n"
			"	z REAL,\n"
			"	FOREIGN KEY (struct_id)\n"
			"		REFERENCES structures (struct_id)\n"
			"		DEFERRABLE INITIALLY DEFERRED);";
	}else if(db_mode == "mysql")
	{
		return
			"CREATE TABLE IF NOT EXISTS nonprotein_residue_conformation (\n"
			"	struct_id INTEGER,\n"
			"	seqpos INTEGER,\n"
			"	phi DOUBLE,\n"
			"	psi DOUBLE,\n"
			"	omega DOUBLE,\n"
			"	FOREIGN KEY (struct_id) REFERENCES structures (struct_id));"
			"\n"
			"CREATE TABLE IF NOT EXISTS nonprotein_residue_angles (\n"
			"	struct_id INTEGER,\n"
			"	seqpos INTEGER,\n"
			"	chinum INTEGER,\n"
			"	chiangle DOUBLE,\n"
			"	FOREIGN KEY (struct_id) REFERENCES structures (struct_id));"
			"\n"
			"CREATE TABLE IF NOT EXISTS residue_atom_coords (\n"
			"	struct_id INTEGER,\n"
			"	seqpos INTEGER,\n"
			"	atomno INTEGER,\n"
			"	x DOUBLE,\n"
			"	y DOUBLE,\n"
			"	z DOUBLE,\n"
			"	FOREIGN KEY (struct_id) REFERENCES structures (struct_id));";
	}else
	{
		return "";
	}

}

Size
ResidueConformationFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size struct_id,
	sessionOP db_session
){
	//check to see if this structure is ideal
	bool ideal = true;
	core::conformation::Conformation const & conformation(pose.conformation());
	for(core::Size resn=1; resn <= pose.n_residue();++resn){
		if(!relevant_residues[resn]) continue;
		bool residue_status(core::conformation::is_ideal_position(resn,conformation));
		if(!residue_status){
			ideal = false;
			break;
		}
	}

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

		statement stmt = (*db_session) <<
			"INSERT INTO nonprotein_residue_conformation VALUES (?,?,?,?,?)" <<
			struct_id << i << phi << psi << omega ;
		stmt.exec();
	
		for(core::Size chi_num = 1; chi_num <= resi.nchi();++chi_num){
			core::Real chi_angle = resi.chi(chi_num);
			while(true)
			{
				try
				{
					statement stmt = (*db_session) <<
							"INSERT INTO nonprotein_residue_angles VALUES (?,?,?,?)" <<
							struct_id << i << chi_num << chi_angle;
					stmt.exec();
					break;
				}catch(cppdb::cppdb_error &)
				{
					usleep(10);
					continue;
				}
			}
		}
		if(!ideal || resi.is_ligand()){ // always store coords for a ligand
			for(Size atom = 1; atom <= resi.natoms(); ++atom){
				core::Vector coords = resi.xyz(atom);
				while(true)
				{
					try
					{
						statement stmt = (*db_session) <<
							"INSERT INTO residue_atom_coords VALUES(?,?,?,?,?,?)" <<
							struct_id << i << atom << coords.x() << coords.y() <<coords.z();
						stmt.exec();
						break;
					}catch(cppdb::cppdb_error &)
					{
						usleep(10);
						continue;
					}
				}
			}
		}
	}
	//transact_guard.commit();
	return 0;
}

void
ResidueConformationFeatures::delete_record(
	Size struct_id,
	sessionOP db_session
){
	while(true)
	{
		try
		{
			statement stmt = (*db_session) << "DELETE FROM nonprotein_residue_conformation WHERE struct_id == ?;\n"<<struct_id;
			stmt.exec();
			stmt = (*db_session) << "DELETE FROM nonprotein_residue_angles WHERE struct_id == ?;\n" << struct_id;
			stmt.exec();
			stmt = (*db_session) << "DELETE FROM residue_atom_coords WHERE struct_id == ?;" << struct_id;
			stmt.exec();
			break;
		}catch(cppdb::cppdb_error &)
		{
			usleep(10);
			continue;
		}
	}
}

void
ResidueConformationFeatures::load_into_pose(
	sessionOP db_session,
	Size struct_id,
	Pose & pose
){
	load_conformation(db_session, struct_id, pose);
}

void
ResidueConformationFeatures::load_conformation(
	sessionOP db_session,
	Size struct_id,
	Pose & pose
){
	if(pose.is_fullatom()){
		while(true)
		{
			try
			{
				result protein_res = (*db_session) <<
					"SELECT\n"
					"	seqpos,\n"
					"	phi,\n"
					"	psi,\n"
					"	omega\n"
					"FROM\n"
					"	nonprotein_residue_conformation\n"
					"WHERE\n"
					"	nonprotein_residue_conformation.struct_id=?;" << struct_id;
				result res_conformation = (*db_session) <<
					"SELECT\n"
					"	seqpos,\n"
					"	chinum,\n"
					"	chiangle\n"
					"FROM\n"
					"	nonprotein_residue_angles\n"
					"WHERE\n"
					"	nonprotein_residue_angles.struct_id=?;" << struct_id;

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
				while(res_conformation.next()){
					//Size nchi(pose.residue_type(seqpos).nchi());
					Size seqpos;
					Size chinum;
					Real chiangle;
					set_coords_for_residue(db_session,struct_id,seqpos,pose);
					res_conformation >> seqpos >> chinum >> chiangle;
					pose.set_chi(chinum,seqpos,chiangle);
				}
				break;
			}catch(cppdb::cppdb_error &)
			{
				usleep(10);
				continue;
			}
		}

	}else{
		while(true)
		{
			try
			{
				result protein_res = (*db_session) <<
					"SELECT\n"
					"	seqpos,\n"
					"	phi,\n"
					"	psi,\n"
					"	omega\n"
					"FROM\n"
					"	nonprotein_residue_conformation\n"
					"WHERE\n"
					"	nonprotein_residue_conformation.struct_id=?;" << struct_id;
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
				break;
			}catch(cppdb::cppdb_error &)
			{
				usleep(10);
				continue;
			}
		}
	}
}

//This should be factored out into a non-member function.
void ResidueConformationFeatures::set_coords_for_residue(
		sessionOP db_session,
		Size struct_id,
		Size seqpos,
		Pose & pose
){

	while(true)
	{
		try
		{
			result res = (*db_session) <<
					"SELECT\n"
					"	atomno,\n"
					"	x,\n"
					"	y,\n"
					"	z\n"
					"FROM\n"
					"	residue_atom_coords\n"
					"WHERE\n"
					"residue_atom_coords.struct_id=? AND residue_atom_coords.seqpos=?;" << struct_id <<	seqpos;
			while(res.next()){
				Size atomno;
				Real x,y,z;
				res >> atomno >> x >> y >> z;

				core::id::AtomID atom_id(atomno,seqpos);
				core::Vector coords(x,y,z);
				pose.set_xyz(atom_id,coords);
			}
			break;
		}catch(cppdb::cppdb_error &)
		{
			usleep(10);
			continue;
		}
	}

}


} // features
} // protocols
