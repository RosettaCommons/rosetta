// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ProteinResidueConformationFeatures.cc
/// @brief  report idealized torsional DOFs Statistics Scientific Benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/ProteinResidueConformationFeatures.hh>

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
#include <utility/string_util.hh>

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
using utility::sql_database::sessionOP;
using utility::vector1;
using cppdb::statement;
using cppdb::result;

string
ProteinResidueConformationFeatures::type_name() const {
	return "ProteinResidueConformationFeatures";
}

string
ProteinResidueConformationFeatures::schema() const {

	std::string db_mode(basic::options::option[basic::options::OptionKeys::inout::database_mode]);
	if(db_mode == "sqlite3")
	{
		return
			"CREATE TABLE IF NOT EXISTS protein_residue_conformation (\n"
			"	struct_id INTEGER,\n"
			"	seqpos INTEGER,\n"
			"	secstruct STRING,\n"
			"	phi REAL,\n"
			"	psi REAL,\n"
			"	omega REAL,\n"
			"	chi1 REAL,\n"
			"	chi2 REAL,\n"
			"	chi3 REAL,\n"
			"	chi4 REAL,\n"
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
			"CREATE TABLE IF NOT EXISTS protein_residue_conformation (\n"
			"	struct_id INTEGER,\n"
			"	seqpos INTEGER,\n"
			"	secstruct TEXT,\n"
			"	phi DOUBLE,\n"
			"	psi DOUBLE,\n"
			"	omega DOUBLE,\n"
			"	chi1 DOUBLE,\n"
			"	chi2 DOUBLE,\n"
			"	chi3 DOUBLE,\n"
			"	chi4 DOUBLE,\n"
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
ProteinResidueConformationFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size struct_id,
	sessionOP db_session
){
	bool fullatom(pose.is_fullatom());

	//check to see if this structure is ideal
	bool ideal = true;
	core::conformation::Conformation const & conformation(pose.conformation());
	for(core::Size resn=1; resn <= pose.n_residue();++resn)
	{
		if(!relevant_residues[resn]) continue;
		bool residue_status(core::conformation::is_ideal_position(resn,conformation));
		if(!residue_status)
		{
			ideal = false;
			break;
		}
	}

	//cppdb::transaction transact_guard(*db_session);
	for (Size i = 1; i <= pose.total_residue(); ++i) {
		if(!relevant_residues[i]) continue;
		Residue const & resi = pose.residue(i);
		if(resi.aa() > num_canonical_aas)
		{
			continue;
		}
		std::string secstruct = utility::to_string<char>(pose.secstruct(i));

		Real phi = 0.0;
		Real psi = 0.0;
		Real omega = 0.0;

		//If you have a non ideal structure, and you store both cartesian coordinates and backbone
		//chi angles and read them into a pose, most of the backbone oxygens will be placed in correctly
		//I currently have absolutely no idea why this is, but this fixes it.
		//It is worth noting that the current implementation of Binary protein silent files does the same thing

		if(ideal)
		{
			 phi = resi.mainchain_torsion(1);
			 psi = resi.mainchain_torsion(2);
			 omega = resi.mainchain_torsion(3);
		}
		Real chi1 = fullatom && resi.nchi() >= 1 ? resi.chi(1) : 0.0;
		Real chi2 = fullatom && resi.nchi() >= 2 ? resi.chi(2) : 0.0;
		Real chi3 = fullatom && resi.nchi() >= 3 ? resi.chi(3) : 0.0;
		Real chi4 = fullatom && resi.nchi() >= 4 ? resi.chi(4) : 0.0;

		statement stmt = (*db_session) <<
			"INSERT INTO protein_residue_conformation VALUES (?,?,?,?,?,?,?,?,?,?)" <<
			struct_id << i << secstruct << phi << psi << omega << chi1 << chi2 << chi3 << chi4 ;
		stmt.exec();

		if(!ideal)
		{
			for(Size atom = 1; atom <= resi.natoms(); ++atom)
			{
				core::Vector coords = resi.xyz(atom);
				statement stmt = (*db_session) <<
					"INSERT INTO residue_atom_coords VALUES(?,?,?,?,?,?)" <<
					struct_id << i << atom << coords.x() << coords.y() <<coords.z();
					//std::cout <<"*" << i << " " << resi.atom_name(atom) << " " << coords.x() << " " << coords.y() << " "<< coords.z() <<std::endl;
				stmt.exec();
			}
		}
	}
	//transact_guard.commit();

	return 0;
}


void
ProteinResidueConformationFeatures::load_into_pose(
	sessionOP db_session,
	Size struct_id,
	Pose & pose
){
	load_conformation(db_session, struct_id, pose);
}

void
ProteinResidueConformationFeatures::load_conformation(
	sessionOP db_session,
	Size struct_id,
	Pose & pose
){
	if(pose.is_fullatom()){
		result res = (*db_session) <<
			"SELECT\n"
			"	seqpos,\n"
			"	secstruct,\n"
			"	phi,\n"
			"	psi,\n"
			"	omega,\n"
			"	chi1,\n"
			"	chi2,\n"
			"	chi3,\n"
			"	chi4\n"
			"FROM\n"
			"	protein_residue_conformation\n"
			"WHERE\n"
			"	protein_residue_conformation.struct_id=?;" << struct_id;
		while(res.next()){
			Size seqpos;
			Real phi,psi,omega,chi1,chi2,chi3,chi4;
			std::string secstruct;
			res >> seqpos >> secstruct >> phi >> psi >> omega >> chi1 >> chi2 >> chi3 >> chi4 ;
			if (!pose.residue_type(seqpos).is_protein()){
				// WARNING why are you storing non-protein in the ProteinSilentReport?
				continue;
			}
			set_coords_for_residue(db_session,struct_id,seqpos,pose);
			if(!(phi < 0.00001 && psi < 0.00001 && omega < 0.00001) )
			{
				pose.set_phi(seqpos,phi);
				pose.set_psi(seqpos,psi);
				pose.set_omega(seqpos,omega);
			}
			pose.set_secstruct(seqpos,static_cast<char>(secstruct[0]));
			Size nchi(pose.residue_type(seqpos).nchi());
			if(1 <= nchi) pose.set_chi(1, seqpos, chi1);
			if(2 <= nchi) pose.set_chi(2, seqpos, chi2);
			if(3 <= nchi) pose.set_chi(3, seqpos, chi3);
			if(4 <= nchi) pose.set_chi(4, seqpos, chi4);
		}

	}else{
		result res = (*db_session) <<
			"SELECT\n"
			"	seqpos,\n"
			"	phi,\n"
			"	psi,\n"
			"	omega\n"
			"FROM\n"
			"	protein_residue_conformation\n"
			"WHERE\n"
			"	protein_residue_conformation.struct_id=?;" << struct_id;
		while(res.next()){
			Size seqpos;
			Real phi,psi,omega;
			res >> seqpos >> phi >> psi >> omega;
			if (!pose.residue_type(seqpos).is_protein()){
				// WARNING why are you storing non-protein in the ProteinSilentReport?
				continue;
			}

			//pose.set_phi(seqpos,phi);
			//pose.set_psi(seqpos,psi);
			//pose.set_omega(seqpos,omega);
			set_coords_for_residue(db_session,struct_id,seqpos,pose);
		}
	}
}


//This should be factored out into a non-member function.
void ProteinResidueConformationFeatures::set_coords_for_residue(
		sessionOP db_session,
		Size struct_id,
		Size seqpos,
		Pose & pose
){

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
	while(res.next())
	{
		Size atomno;
		Real x,y,z;
		res >> atomno >> x >> y >> z;

		core::id::AtomID atom_id(atomno,seqpos);
		core::Vector coords(x,y,z);
		pose.set_xyz(atom_id,coords);
		//std::cout <<seqpos << " " << pose.residue(seqpos).atom_name(atomno)<< " " << x << " "<<y << " "<<z <<std::endl;
	}

}

} // namespace
} // namespace
