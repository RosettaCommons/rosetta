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
#include <core/chemical/orbitals/OrbitalType.hh>
#include <core/chemical/AtomType.hh>
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

#include <utility/exit.hh>

#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>

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

OrbitalsFeatures::OrbitalsFeatures()
{
	if(!basic::options::option[ basic::options::OptionKeys::in::add_orbitals]){
		utility_exit_with_message( "Trying to run features test without orbitals! Pass the flag -add_orbitals!" );
	}
}

OrbitalsFeatures::OrbitalsFeatures( OrbitalsFeatures const & ) :
			FeaturesReporter()
{
	if(!basic::options::option[ basic::options::OptionKeys::in::add_orbitals]){
		utility_exit_with_message( "Trying to run features test without orbitals! Pass the flag -add_orbitals!" );
	}
}

OrbitalsFeatures::~OrbitalsFeatures(){}

string
OrbitalsFeatures::type_name() const { return "OrbitalsFeatures"; }

string
OrbitalsFeatures::schema() const {
	return
			"CREATE TABLE IF NOT EXISTS HPOL_orbital (\n"
			"	struct_id TEXT,\n"
			"	resNum1 INTEGER,\n"
			"	orbNum1 INTEGER,\n"
			"	orbName1 TEXT,\n"
			"	resNum2 INTEGER,\n"
			"	hpolNum2 INTEGER,\n"
			"   htype2 TEXT,\n"
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
			"CREATE TABLE IF NOT EXISTS HARO_orbital (\n"
			"	struct_id TEXT,\n"
			"	resNum1 INTEGER,\n"
			"	orbNum1 INTEGER,\n"
			"	orbName1 TEXT,\n"
			"	resNum2 INTEGER,\n"
			"	haroNum2 INTEGER,\n"
			"   htype2 TEXT,\n"
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

utility::vector1<std::string>
OrbitalsFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	return dependencies;
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
		vector1< bool > const & relevant_residues,
		Size const struct_id,
		sessionOP db_session
){
	std::string orbita_H_string = "INSERT INTO HPOL_orbital VALUES (?,?,?,?,?,?,?,?,?,?);";
	statement orbital_H_statement(basic::database::safely_prepare_statement(orbita_H_string,db_session));

	for(Size resNum1 = 1; resNum1 <= pose.n_residue(); ++resNum1){
		Residue const res1 = pose.residue(resNum1);
		Size resNum2(0);
		string orbName1;
		string htype2;
		Size orbNum1(0);
		Size hpolNum2(0);
		Size haroNum2(0);
		Size atm(0);
		Real AOH_angle(0);
		Real DHO_angle(0);
		Real dist(10.1); // min distance used to derive statistics. Should be the shortest distance between an orbital and hydrogen
		for(Size res_num2 = 1; res_num2 <= pose.n_residue(); ++res_num2){
			foreach(Size const atm1, res1.atoms_with_orb_index()){
				foreach(Size const temp_orbNum1, res1.bonded_orbitals(atm1)){
					xyzVector<Real> const res1_orb_xyz(res1.orbital_xyz(temp_orbNum1));
					Residue const res2 = pose.residue(res_num2);
					if(resNum1 != res_num2){
						foreach(Size const temp_hpolNum2, res2.Hpol_index()){
							if(res1.atom_is_backbone(atm1) && res2.atom_is_backbone(temp_hpolNum2)){
								continue;//do nothing. Dont really want to calculate bb-bb interactions
							}else{
								xyzVector<Real> res2_H_xyz(res2.atom(temp_hpolNum2).xyz());
								Real const container(res1_orb_xyz.distance(res2_H_xyz));
								if(container <= dist){
									xyzVector<Real> const bonded_atom_to_orbital_xyz(res1.atom(atm1).xyz());
									Size donor_atom(res2.bonded_neighbor(temp_hpolNum2)[1]);
									xyzVector<Real> donor(res2.xyz(donor_atom));
									AOH_angle = cos_of(bonded_atom_to_orbital_xyz, res1_orb_xyz, res2_H_xyz );
									DHO_angle = cos_of(donor, res2_H_xyz, res1_orb_xyz);
									dist=container;
									resNum2=res_num2;
									hpolNum2=temp_hpolNum2;
									htype2=res2.atom_type(temp_hpolNum2).name();
									orbNum1=temp_orbNum1;
									orbName1=res1.orbital_type(temp_orbNum1).name();
									atm=atm1;
								}
							}
						}
					}
				}
			}

		}
		if(dist <=10.0){
			Residue const res2 = pose.residue(resNum2);
			orbital_H_statement.bind(1, struct_id);
			orbital_H_statement.bind(2, resNum1);
			orbital_H_statement.bind(3,orbNum1 );
			orbital_H_statement.bind(4, orbName1);
			orbital_H_statement.bind(5, resNum2);
			orbital_H_statement.bind(6, hpolNum2);
			orbital_H_statement.bind(7, htype2);
			orbital_H_statement.bind(8, dist);
			orbital_H_statement.bind(9, AOH_angle);
			orbital_H_statement.bind(10, DHO_angle);
			basic::database::safely_write_to_database(orbital_H_statement);
		}
	}
}


///@brief get statistics based upon aromatic hydrogen to orbital distance/angle
void
OrbitalsFeatures::report_haro_orbital_interactions(
		Pose const & pose,
		vector1< bool > const & relevant_residues,
		Size const struct_id,
		sessionOP db_session
){
	std::string orbita_H_string = "INSERT INTO HARO_orbital VALUES (?,?,?,?,?,?,?,?,?,?);";
	statement orbital_H_statement(basic::database::safely_prepare_statement(orbita_H_string,db_session));


	for(Size resNum1 = 1; resNum1 <= pose.n_residue(); ++resNum1){
		Residue const res1 = pose.residue(resNum1);
		Size resNum2(0);
		string orbName1;
		string htype2;
		Size orbNum1(0);
		Size haroNum2(0);
		Size atm(0);
		Real AOH_angle(0);
		Real DHO_angle(0);
		Real dist(10.1); // min distance used to derive statistics. Should be the shortest distance between an orbital and hydrogen
		for(Size res_num2 = 1; res_num2 <= pose.n_residue(); ++res_num2){
			foreach(Size const atm1, res1.atoms_with_orb_index()){
				foreach(Size const temp_orbNum1, res1.bonded_orbitals(atm1)){
					xyzVector<Real> const res1_orb_xyz(res1.orbital_xyz(temp_orbNum1));
					Residue const res2 = pose.residue(res_num2);
					if(resNum1 != res_num2){
						foreach(Size const temp_haroNum2, res2.Haro_index()){
							xyzVector<Real> res2_H_xyz(res2.atom(temp_haroNum2).xyz());
							Real const container(res1_orb_xyz.distance(res2_H_xyz));
							if(container <= dist){
								xyzVector<Real> const bonded_atom_to_orbital_xyz(res1.atom(atm1).xyz());
								Size donor_atom(res2.bonded_neighbor(temp_haroNum2)[1]);
								xyzVector<Real> donor(res2.xyz(donor_atom));
								AOH_angle = cos_of(bonded_atom_to_orbital_xyz, res1_orb_xyz, res2_H_xyz );
								DHO_angle = cos_of(donor, res2_H_xyz, res1_orb_xyz);
								dist=container;
								resNum2=res_num2;
								haroNum2=temp_haroNum2;
								htype2=res2.atom_type(temp_haroNum2).name();
								orbNum1=temp_orbNum1;
								orbName1= res1.orbital_type(temp_orbNum1).name();
								atm=atm1;
							}
						}
					}
				}

			}
		}
		if(dist <=10.0){
			orbital_H_statement.bind(1, struct_id);
			orbital_H_statement.bind(2, resNum1);
			orbital_H_statement.bind(3,orbNum1 );
			orbital_H_statement.bind(4, orbName1);
			orbital_H_statement.bind(5, resNum2);
			orbital_H_statement.bind(6, haroNum2);
			orbital_H_statement.bind(7, htype2);
			orbital_H_statement.bind(8, dist);
			orbital_H_statement.bind(9, AOH_angle);
			orbital_H_statement.bind(10, DHO_angle);
			basic::database::safely_write_to_database(orbital_H_statement);
		}
	}
}






} // namesapce
} // namespace

