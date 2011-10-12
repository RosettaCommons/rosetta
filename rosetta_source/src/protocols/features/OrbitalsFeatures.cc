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
		"CREATE TABLE IF NOT EXISTS orbital_polar_hydrogen_interactions (\n"
		"	struct_id TEXT,\n"
		"	resNum1 INTEGER,\n"
		"	orbNum1 INTEGER,\n"
		"	orbName1 TEXT,\n"
		"	resNum2 INTEGER,\n"
		"	hpolNum2 INTEGER,\n"
		"	dist REAL,\n"
		"	angle REAL,\n"
		"	FOREIGN KEY (struct_id, resNum1)\n"
		"		REFERENCES residues (struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	FOREIGN KEY (struct_id, resNum2)\n"
		"		REFERENCES residues (struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	PRIMARY KEY(struct_id, resNum1, orbNum1, resNum2, hpolNum2));\n"
		"\n"
		"CREATE TABLE IF NOT EXISTS orbital_aromatic_hydrogen_interactions (\n"
		"	struct_id TEXT,\n"
		"	resNum1 INTEGER,\n"
		"	orbNum1 INTEGER,\n"
		"	orbName1 TEXT,\n"
		"	resNum2 INTEGER,\n"
		"	haroNum2 INTEGER,\n"
		"	dist REAL,\n"
		"	angle REAL,\n"
		"	FOREIGN KEY (struct_id, resNum1)\n"
		"		REFERENCES residues (struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	FOREIGN KEY (struct_id, resNum2)\n"
		"		REFERENCES residues (struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	PRIMARY KEY(struct_id, resNum1, orbNum1, resNum2, haroNum2));\n"
		"\n"
		"CREATE TABLE IF NOT EXISTS orbital_orbital_interactions (\n"
		"	struct_id TEXT,\n"
		"	resNum1 INTEGER,\n"
		"	orbNum1 INTEGER,\n"
		"	orbName1 TEXT,\n"
		"	resNum2 INTEGER,\n"
		"	orbNum2 INTEGER,\n"
		"	dist REAL,\n"
		"	angle REAL,\n"
		"	FOREIGN KEY (struct_id, resNum1)\n"
		"		REFERENCES residues (struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	FOREIGN KEY (struct_id, resNum2)\n"
		"		REFERENCES residues (struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	PRIMARY KEY(struct_id, resNum1, orbNum1, resNum2, orbNum2));\n";


}

Size
OrbitalsFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size const struct_id,
	sessionOP db_session
){
	report_orbital_interactions( pose, relevant_residues, struct_id, db_session );
	return 0;
}


///@brief get statistics based upon hydrogen to orbital distance/angle
void
OrbitalsFeatures::report_orbital_interactions(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size const struct_id,
	sessionOP db_session
){
	// This should match OrbitalsScore::atomic_interaction_cutoff()
	Real atomic_interaction_cutoff = 4.5; // min distance

	for(Size resNum1 = 1; resNum1 <= pose.n_residue(); ++resNum1){
		if(!relevant_residues[resNum1]) continue;
		Residue const res1 = pose.residue(resNum1);
		foreach(Size const atm1, res1.atoms_with_orb_index()){
			xyzVector<Real> const bonded_atom_xyz(res1.atom(atm1).xyz());
			foreach(Size const orbNum1, res1.bonded_orbitals(atm1)){
				string const orbName1(res1.mm_atom_name(orbNum1));
				xyzVector<Real> const res1_orb_xyz(res1.orbital_xyz(orbNum1));

				for(Size resNum2 = 1; resNum2 <= pose.n_residue(); ++resNum2){
					if(!relevant_residues[resNum2]) continue;
					Residue const res2 = pose.residue(resNum2);
					if(resNum1 != resNum2){
						foreach(Size const hpolNum2, res2.Hpos_polar_sc()){
							xyzVector<Real> res2_H_xyz(res2.atom(hpolNum2).xyz());
							Real const dist(res1_orb_xyz.distance(res2_H_xyz));
							if(dist <= atomic_interaction_cutoff){
								Real angle = cos_of(bonded_atom_xyz, res1_orb_xyz, res2_H_xyz );
								statement stmt = (*db_session)
									<< "INSERT INTO orbital_polar_hydrogen_interactions VALUES (?,?,?,?,?,?,?,?);"
									<< struct_id
									<< resNum1 << orbNum1 << orbName1
									<< resNum2 << hpolNum2
									<< dist << angle;
								basic::database::safely_write_to_database(stmt);
							}
						}
						foreach(Size const haroNum2, res2.Haro_index()){
							xyzVector<Real> res2_H_xyz(res2.atom(haroNum2).xyz());
							Real const dist(res1_orb_xyz.distance(res2_H_xyz));
							if(dist <= atomic_interaction_cutoff){
								Real angle = cos_of(bonded_atom_xyz, res1_orb_xyz, res2_H_xyz );
								statement stmt = (*db_session)
									<< "INSERT INTO orbital_aromatic_hydrogen_interactions VALUES (?,?,?,?,?,?,?,?);"
									<< struct_id
									<< resNum1 << orbNum1 << orbName1
									<< resNum2 << haroNum2
									<< dist << angle;
								basic::database::safely_write_to_database(stmt);
							}
						}
						if(resNum1 < resNum2){
							for(Size orbNum2=1; orbNum2 <= res2.n_orbitals(); ++orbNum2){
								xyzVector<Real> res2_orb_xyz(res2.orbital_xyz(orbNum2));
								Real const dist(res1_orb_xyz.distance(res2_orb_xyz));
								if(dist <= atomic_interaction_cutoff){
									Real angle = cos_of(bonded_atom_xyz, res1_orb_xyz, res2_orb_xyz );
									statement stmt = (*db_session)
										<< "INSERT INTO orbital_orbital_interactions VALUES (?,?,?,?,?,?,?,?);"
										<< struct_id
										<< resNum1 << orbNum1 << orbName1
										<< resNum2 << orbNum2
										<< dist << angle;
									basic::database::safely_write_to_database(stmt);

								}
							}
						}
					}
				}
			}
		}
	}
}



} // namesapce
} // namespace
