// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/OrbitalsFeatures.cc
/// @brief  report orbital geometry and scores to features statistics scientific benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/OrbitalsFeatures.hh>

// Project Headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/orbitals/OrbitalsLookup.hh>
#include <core/types.hh>
#include <basic/database/sql_utils.hh>

//Numeric Headers
#include <numeric/xyzVector.hh>

//Utility Headers
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

// Boost Headers
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// External Headers
#include <cppdb/frontend.h>

namespace protocols{
namespace features{

using std::string;
using core::chemical::AtomIndices;
using core::pose::Pose;
using core::Size;
using core::Real;
using core::conformation::Residue;
using core::scoring::orbitals::OrbitalsLookup;
using numeric::xyzVector;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;

OrbitalsFeatures::OrbitalsFeatures(){}

OrbitalsFeatures::OrbitalsFeatures( OrbitalsFeatures const & ) :
	FeaturesReporter()
{}

OrbitalsFeatures::~OrbitalsFeatures(){}

string
OrbitalsFeatures::type_name() const { return "OrbitalsFeatures"; }

string
OrbitalsFeatures::schema() const {
	return
		"CREATE TABLE IF NOT EXISTS HPOL_sc_H_sc_orb (\n"
		"	struct_id TEXT,\n"
		"	resNum1 INTEGER,\n"
		"	orbNum1 INTEGER,\n"
		"	orbName1 TEXT,\n"
		"	resNum2 INTEGER,\n"
		"	hpolNum2 INTEGER,\n"
		"	dist REAL,\n"
		"	AOH_angle REAL,\n"
		"   DHO_angle REAL, \n"
		"	FOREIGN KEY (struct_id, resNum1)\n"
		"		REFERENCES residues (struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	FOREIGN KEY (struct_id, resNum2)\n"
		"		REFERENCES residues (struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	PRIMARY KEY(struct_id, resNum1, orbName1, resNum2, hpolNum2));\n"
		"\n"
		"CREATE TABLE IF NOT EXISTS HPOL_bb_H_sc_orb (\n"
		"	struct_id TEXT,\n"
		"	resNum1 INTEGER,\n"
		"	orbNum1 INTEGER,\n"
		"	orbName1 TEXT,\n"
		"	resNum2 INTEGER,\n"
		"	hpolNum2 INTEGER,\n"
		"	dist REAL,\n"
		"	AOH_angle REAL,\n"
		"   DHO_angle REAL, \n"
		"	FOREIGN KEY (struct_id, resNum1)\n"
		"		REFERENCES residues (struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	FOREIGN KEY (struct_id, resNum2)\n"
		"		REFERENCES residues (struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	PRIMARY KEY(struct_id, resNum1, orbName1, resNum2, hpolNum2));\n"
		"\n"
		"CREATE TABLE IF NOT EXISTS HPOL_sc_H_bb_orb (\n"
		"	struct_id TEXT,\n"
		"	resNum1 INTEGER,\n"
		"	orbNum1 INTEGER,\n"
		"	orbName1 TEXT,\n"
		"	resNum2 INTEGER,\n"
		"	hpolNum2 INTEGER,\n"
		"	dist REAL,\n"
		"	AOH_angle REAL,\n"
		"   DHO_angle REAL, \n"
		"	FOREIGN KEY (struct_id, resNum1)\n"
		"		REFERENCES residues (struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	FOREIGN KEY (struct_id, resNum2)\n"
		"		REFERENCES residues (struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	PRIMARY KEY(struct_id, resNum1, orbName1, resNum2, hpolNum2));\n"
		"\n"
		"CREATE TABLE IF NOT EXISTS HARO_sc_H_sc_orb (\n"
		"	struct_id TEXT,\n"
		"	resNum1 INTEGER,\n"
		"	orbNum1 INTEGER,\n"
		"	orbName1 TEXT,\n"
		"	resNum2 INTEGER,\n"
		"	haroNum2 INTEGER,\n"
		"	dist REAL,\n"
		"	AOH_angle REAL,\n"
		"   DHO_angle REAL, \n"
		"	FOREIGN KEY (struct_id, resNum1)\n"
		"		REFERENCES residues (struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	FOREIGN KEY (struct_id, resNum2)\n"
		"		REFERENCES residues (struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	PRIMARY KEY(struct_id, resNum1, orbName1, resNum2, haroNum2));\n"
		"\n";


}

Size
OrbitalsFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size const struct_id,
	sessionOP db_session
){
	report_hpol_orbital_interactions( pose, relevant_residues, struct_id, db_session );
	report_haro_orbital_interactions( pose, relevant_residues, struct_id, db_session );
	return 0;
}


///@brief get statistics based upon polar hydrogen to orbital distance/angle
void
OrbitalsFeatures::report_hpol_orbital_interactions(
	Pose const & pose,
	vector1< bool > const & /*relevant_residues*/,
	Size const struct_id,
	sessionOP db_session
){
	//should match max_dist_squared_^1/2

	std::string sc_bb_string = "INSERT INTO HPOL_sc_H_bb_orb VALUES (?,?,?,?,?,?,?,?,?);";
	std::string bb_sc_string = "INSERT INTO HPOL_bb_H_sc_orb VALUES (?,?,?,?,?,?,?,?,?);";
	std::string sc_sc_string = "INSERT INTO HPOL_sc_H_sc_orb VALUES (?,?,?,?,?,?,?,?,?);";

	statement sc_bb_statement(basic::database::safely_prepare_statement(sc_bb_string,db_session));
	statement bb_sc_statement(basic::database::safely_prepare_statement(bb_sc_string,db_session));
	statement sc_sc_statement(basic::database::safely_prepare_statement(sc_sc_string,db_session));

	for(Size resNum1 = 1; resNum1 <= pose.n_residue(); ++resNum1){
		Residue const res1 = pose.residue(resNum1);
		foreach(Size const atm1, res1.atoms_with_orb_index()){
			Size resNum2(0);
			string orbName1;
			Size orbNum1(0);
			Size hpolNum2(0);
			Size haroNum2(0);
			Size atm(0);
			Real AOH_angle(0);
			Real DHO_angle(0);
			Real dist(10.1); // min distance used to derive statistics. Should be the shortest distance between an orbital and hydrogen
			xyzVector<Real> const bonded_atom_to_orbital_xyz(res1.atom(atm1).xyz());
			foreach(Size const temp_orbNum1, res1.bonded_orbitals(atm1)){
				xyzVector<Real> const res1_orb_xyz(res1.orbital_xyz(temp_orbNum1));
				for(Size res_num2 = 1; res_num2 <= pose.n_residue(); ++res_num2){
					Residue const res2 = pose.residue(res_num2);
					if(resNum1 != res_num2){
						foreach(Size const temp_hpolNum2, res2.Hpol_index()){
							if(res1.atom_is_backbone(atm1) && res2.atom_is_backbone(temp_hpolNum2)){
								continue;//do nothing. Dont really want to calculate bb-bb interactions
							}else{
								xyzVector<Real> res2_H_xyz(res2.atom(temp_hpolNum2).xyz());
								Real const container(res1_orb_xyz.distance(res2_H_xyz));
								if(container <= dist){
									Size donor_atom(res2.bonded_neighbor(temp_hpolNum2)[1]);
									xyzVector<Real> donor(res2.xyz(donor_atom));
									AOH_angle = cos_of(bonded_atom_to_orbital_xyz, res1_orb_xyz, res2_H_xyz );
									DHO_angle = cos_of(donor, res2_H_xyz, res1_orb_xyz);
									dist=container;
									resNum2=res_num2;
									hpolNum2=temp_hpolNum2;
									orbNum1=temp_orbNum1;
									orbName1=res1.orbital_name(temp_orbNum1);
									atm=atm1;
								}
							}
						}
					}
				}
				}
			if(dist <=10.0){
				Residue const res2 = pose.residue(resNum2);
					if(res1.atom_is_backbone(atm) && !res2.atom_is_backbone(hpolNum2)){//HPOL_sc_H_bb_orb

						sc_bb_statement.bind(1,struct_id);
						sc_bb_statement.bind(2,resNum1);
						sc_bb_statement.bind(3,orbNum1);
						sc_bb_statement.bind(4,orbName1);
						sc_bb_statement.bind(5,resNum2);
						sc_bb_statement.bind(6,hpolNum2);
						sc_bb_statement.bind(7,dist);
						sc_bb_statement.bind(8,AOH_angle);
						sc_bb_statement.bind(9,DHO_angle);
						basic::database::safely_write_to_database(sc_bb_statement);
					}
					else if(!res1.atom_is_backbone(atm) && res2.atom_is_backbone(hpolNum2)){//HPOL_bb_H_sc_orb

						bb_sc_statement.bind(1,struct_id);
						bb_sc_statement.bind(2,resNum1);
						bb_sc_statement.bind(3,orbNum1);
						bb_sc_statement.bind(4,orbName1);
						bb_sc_statement.bind(5,resNum2);
						bb_sc_statement.bind(6,hpolNum2);
						bb_sc_statement.bind(7,dist);
						bb_sc_statement.bind(8,AOH_angle);
						bb_sc_statement.bind(9,DHO_angle);
						basic::database::safely_write_to_database(bb_sc_statement);
					}
					else if(!res1.atom_is_backbone(atm) && !res2.atom_is_backbone(hpolNum2)){//HPOL_sc_H_sc_orb

						sc_sc_statement.bind(1,struct_id);
						sc_sc_statement.bind(2,resNum1);
						sc_sc_statement.bind(3,orbNum1);
						sc_sc_statement.bind(4,orbName1);
						sc_sc_statement.bind(5,resNum2);
						sc_sc_statement.bind(6,hpolNum2);
						sc_sc_statement.bind(7,dist);
						sc_sc_statement.bind(8,AOH_angle);
						sc_sc_statement.bind(9,DHO_angle);
						basic::database::safely_write_to_database(sc_sc_statement);
					}

			}
		}
	}
}


	///@brief get statistics based upon aromatic hydrogen to orbital distance/angle
void
OrbitalsFeatures::report_haro_orbital_interactions(
	Pose const & pose,
	vector1< bool > const & /*relevant_residues*/,
	Size const struct_id,
	sessionOP db_session
){
	std::string statement_string = "INSERT INTO HARO_sc_H_sc_orb VALUES (?,?,?,?,?,?,?,?,?);";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	for(Size resNum1 = 1; resNum1 <= pose.n_residue(); ++resNum1){
		Residue const res1 = pose.residue(resNum1);
		foreach(Size const atm1, res1.atoms_with_orb_index()){
			Size resNum2(0);
			string orbName1;
			Size orbNum1(0);
			Size haroNum2(0);
			Size atm(0);
			Real AOH_angle(0);
			Real DHO_angle(0);
			Real dist(10.1); // min distance used to derive statistics. Should be the shortest distance between an orbital and hydrogen
			xyzVector<Real> const bonded_atom_to_orbital_xyz(res1.atom(atm1).xyz());
			foreach(Size const temp_orbNum1, res1.bonded_orbitals(atm1)){
				xyzVector<Real> const res1_orb_xyz(res1.orbital_xyz(temp_orbNum1));
				for(Size res_num2 = 1; res_num2 <= pose.n_residue(); ++res_num2){
					Residue const res2 = pose.residue(res_num2);
					if(resNum1 != res_num2){
						foreach(Size const temp_haroNum2, res2.Haro_index()){
							xyzVector<Real> res2_H_xyz(res2.atom(temp_haroNum2).xyz());
							Real const container(res1_orb_xyz.distance(res2_H_xyz));
							if(container <= dist){
								Size donor_atom(res2.bonded_neighbor(temp_haroNum2)[1]);
								xyzVector<Real> donor(res2.xyz(donor_atom));
								AOH_angle = cos_of(bonded_atom_to_orbital_xyz, res1_orb_xyz, res2_H_xyz );
								DHO_angle = cos_of(donor, res2_H_xyz, res1_orb_xyz);
								dist=container;
								resNum2=res_num2;
								haroNum2=temp_haroNum2;
								orbNum1=temp_orbNum1;
								orbName1=res1.orbital_name(temp_orbNum1);
								atm=atm1;
							}
						}
					}
				}

			}
			if(dist <=10.0){
				Residue const res2 = pose.residue(resNum2);

				stmt.bind(1,struct_id);
				stmt.bind(2,resNum1);
				stmt.bind(3,orbNum1);
				stmt.bind(4,orbName1);
				stmt.bind(5,resNum2);
				stmt.bind(6,haroNum2);
				stmt.bind(7,dist);
				stmt.bind(8,AOH_angle);
				stmt.bind(9,DHO_angle);
				basic::database::safely_write_to_database(stmt);
			}
		}
	}
}






} // namesapce
} // namespace
