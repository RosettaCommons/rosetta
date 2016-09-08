// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>
#include <core/chemical/orbitals/OrbitalType.hh>
#include <core/chemical/AtomType.hh>
//Numeric Headers
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>


//Utility Headers
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

// Boost Headers
#include <boost/foreach.hpp>

// External Headers
#include <cppdb/frontend.h>

#include <utility/exit.hh>

#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>

namespace protocols {
namespace features {

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
	if ( basic::options::option[ basic::options::OptionKeys::in::add_orbitals] != 1 ) {
		utility_exit_with_message( "Trying to run features test without orbitals! Pass the flag -add_orbitals!" );
	}
}

OrbitalsFeatures::OrbitalsFeatures( OrbitalsFeatures const & ) :
	FeaturesReporter()
{
	if ( basic::options::option[ basic::options::OptionKeys::in::add_orbitals] != 1 ) {
		utility_exit_with_message( "Trying to run features test without orbitals! Pass the flag -add_orbitals!" );
	}
}

OrbitalsFeatures::~OrbitalsFeatures()= default;

string
OrbitalsFeatures::type_name() const { return "OrbitalsFeatures"; }
void
OrbitalsFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	write_HPOL_orbital_table_schema(db_session);
	write_HARO_orbital_table_schema(db_session);
	write_orbital_orbital_table_schema(db_session);
}

void
OrbitalsFeatures::write_HPOL_orbital_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column resNum1("resNum1", DbDataTypeOP( new DbInteger() ));
	Column resName1("resName1", DbDataTypeOP( new DbText() ));
	Column orbNum1("orbNum1", DbDataTypeOP( new DbInteger() ));
	Column orbName1("orbName1", DbDataTypeOP( new DbText() ));
	Column resNum2("resNum2", DbDataTypeOP( new DbInteger() ));
	Column resName2("resName2", DbDataTypeOP( new DbText() ));
	Column hpolNum2("hpolNum2", DbDataTypeOP( new DbInteger() ));
	Column htype2("htype2", DbDataTypeOP( new DbText() ));
	Column OrbHdist("OrbHdist", DbDataTypeOP( new DbReal() ));
	Column cosAOH("cosAOH", DbDataTypeOP( new DbReal() ));
	Column cosDHO("cosDHO", DbDataTypeOP( new DbReal() ));
	Column chiBAOH("chiBAOH", DbDataTypeOP( new DbReal() ));
	Column chiBDHO("chiBDHO", DbDataTypeOP( new DbReal() ));
	Column AOH_angle("AOH_angle", DbDataTypeOP( new DbReal() ));
	Column DHO_angle("DHO_angle", DbDataTypeOP( new DbReal() ));
	Column chiBAHD("chiBAHD", DbDataTypeOP( new DbReal() ));
	Column cosAHD("cosAHD", DbDataTypeOP( new DbReal() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(resNum1);
	primary_key_columns.push_back(orbName1);
	primary_key_columns.push_back(resNum2);
	primary_key_columns.push_back(hpolNum2);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns1;
	foreign_key_columns1.push_back(struct_id);
	foreign_key_columns1.push_back(resNum1);
	vector1< std::string > reference_columns1;
	reference_columns1.push_back("struct_id");
	reference_columns1.push_back("resNum");
	ForeignKey foreign_key1(foreign_key_columns1, "residues", reference_columns1, true);

	Columns foreign_key_columns2;
	foreign_key_columns2.push_back(struct_id);
	foreign_key_columns2.push_back(resNum2);
	vector1< std::string > reference_columns2;
	reference_columns2.push_back("struct_id");
	reference_columns2.push_back("resNum");
	ForeignKey foreign_key2(foreign_key_columns2, "residues", reference_columns2, true);


	Schema table("HPOL_orbital", primary_key);
	table.add_foreign_key(foreign_key1);
	table.add_foreign_key(foreign_key2);
	table.add_column(resName1);
	table.add_column(orbNum1);
	table.add_column(resName2);
	table.add_column(hpolNum2);
	table.add_column(htype2);
	table.add_column(OrbHdist);
	table.add_column(cosAOH);
	table.add_column(cosDHO);
	table.add_column(chiBAOH);
	table.add_column(chiBDHO);
	table.add_column(AOH_angle);
	table.add_column(DHO_angle);
	table.add_column(chiBAHD);
	table.add_column(cosAHD);

	table.write(db_session);
}

void
OrbitalsFeatures::write_HARO_orbital_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column resNum1("resNum1", DbDataTypeOP( new DbInteger() ));
	Column resName1("resName1", DbDataTypeOP( new DbText() ));
	Column orbNum1("orbNum1", DbDataTypeOP( new DbInteger() ));
	Column orbName1("orbName1", DbDataTypeOP( new DbText() ));
	Column resNum2("resNum2", DbDataTypeOP( new DbInteger() ));
	Column resName2("resName2", DbDataTypeOP( new DbText() ));
	Column haroNum2("haroNum2", DbDataTypeOP( new DbInteger() ));
	Column htype2("htype2", DbDataTypeOP( new DbText() ));
	Column OrbHdist("orbHdist", DbDataTypeOP( new DbReal() ));
	Column cosAOH("cosAOH", DbDataTypeOP( new DbReal() ));
	Column cosDHO("cosDHO", DbDataTypeOP( new DbReal() ));
	Column chiBAOH("chiBAOH", DbDataTypeOP( new DbReal() ));
	Column chiBDHO("chiBDHO", DbDataTypeOP( new DbReal() ));
	Column AOH_angle("AOH_angle", DbDataTypeOP( new DbReal() ));
	Column DHO_angle("DHO_angle", DbDataTypeOP( new DbReal() ));
	Column chiBAHD("chiBAHD", DbDataTypeOP( new DbReal() ));
	Column cosAHD("cosAHD", DbDataTypeOP( new DbReal() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(resNum1);
	primary_key_columns.push_back(orbName1);
	primary_key_columns.push_back(resNum2);
	primary_key_columns.push_back(haroNum2);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns1;
	foreign_key_columns1.push_back(struct_id);
	foreign_key_columns1.push_back(resNum1);
	vector1< std::string > reference_columns1;
	reference_columns1.push_back("struct_id");
	reference_columns1.push_back("resNum");
	ForeignKey foreign_key1(foreign_key_columns1, "residues", reference_columns1, true);

	Columns foreign_key_columns2;
	foreign_key_columns2.push_back(struct_id);
	foreign_key_columns2.push_back(resNum2);
	vector1< std::string > reference_columns2;
	reference_columns2.push_back("struct_id");
	reference_columns2.push_back("resNum");
	ForeignKey foreign_key2(foreign_key_columns2, "residues", reference_columns2, true);


	Schema table("HARO_orbital", primary_key);
	table.add_foreign_key(foreign_key1);
	table.add_foreign_key(foreign_key2);
	table.add_column(resName1);
	table.add_column(orbNum1);
	table.add_column(resName2);
	table.add_column(haroNum2);
	table.add_column(htype2);
	table.add_column(OrbHdist);
	table.add_column(cosAOH);
	table.add_column(cosDHO);
	table.add_column(chiBAOH);
	table.add_column(chiBDHO);
	table.add_column(AOH_angle);
	table.add_column(DHO_angle);
	table.add_column(chiBAHD);
	table.add_column(cosAHD);

	table.write(db_session);
}

void
OrbitalsFeatures::write_orbital_orbital_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column resNum1("resNum1", DbDataTypeOP( new DbInteger() ));
	Column resName1("resName1", DbDataTypeOP( new DbText() ));
	Column orbNum1("orbNum1", DbDataTypeOP( new DbInteger() ));
	Column orbName1("orbName1", DbDataTypeOP( new DbText() ));
	Column resNum2("resNum2", DbDataTypeOP( new DbInteger() ));
	Column resName2("resName2", DbDataTypeOP( new DbText() ));
	Column orbNum2("orbNum2", DbDataTypeOP( new DbInteger() ));
	Column orbName2("orbName2", DbDataTypeOP( new DbText() ));
	Column orbOrbdist("orbOrbdist", DbDataTypeOP( new DbReal() ));
	Column cosAOO("cosAOO", DbDataTypeOP( new DbReal() ));
	Column cosDOO("cosDOO", DbDataTypeOP( new DbReal() ));
	Column chiBAOO("chiBAOO", DbDataTypeOP( new DbReal() ));
	Column chiBDOO("chiBDOO", DbDataTypeOP( new DbReal() ));
	Column AOO_angle("AOO_angle", DbDataTypeOP( new DbReal() ));
	Column DOO_angle("DOO_angle", DbDataTypeOP( new DbReal() ));
	Column DOA_angle("DOA_angle", DbDataTypeOP( new DbReal() ));
	Column AOD_angle("AOD_angle", DbDataTypeOP( new DbReal() ));
	Column chiBAHD("chiBAHD", DbDataTypeOP( new DbReal() ));
	Column cosAHD("cosAHD", DbDataTypeOP( new DbReal() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(resNum1);
	primary_key_columns.push_back(orbName1);
	primary_key_columns.push_back(resNum2);
	primary_key_columns.push_back(orbNum2);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns1;
	foreign_key_columns1.push_back(struct_id);
	foreign_key_columns1.push_back(resNum1);
	vector1< std::string > reference_columns1;
	reference_columns1.push_back("struct_id");
	reference_columns1.push_back("resNum");
	ForeignKey foreign_key1(foreign_key_columns1, "residues", reference_columns1, true);

	Columns foreign_key_columns2;
	foreign_key_columns2.push_back(struct_id);
	foreign_key_columns2.push_back(resNum2);
	vector1< std::string > reference_columns2;
	reference_columns2.push_back("struct_id");
	reference_columns2.push_back("resNum");
	ForeignKey foreign_key2(foreign_key_columns2, "residues", reference_columns2, true);


	Schema table("orbital_orbital", primary_key);
	table.add_foreign_key(foreign_key1);
	table.add_foreign_key(foreign_key2);
	table.add_column(resName1);
	table.add_column(orbNum1);
	table.add_column(resName2);
	table.add_column(orbName2);
	table.add_column(orbOrbdist);
	table.add_column(cosAOO);
	table.add_column(cosDOO);
	table.add_column(chiBAOO);
	table.add_column(chiBDOO);
	table.add_column(AOO_angle);
	table.add_column(DOO_angle);
	table.add_column(DOA_angle);
	table.add_column(AOD_angle);
	table.add_column(chiBAHD);
	table.add_column(cosAHD);

	table.write(db_session);
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
	StructureID const struct_id,
	sessionOP db_session
){
	report_hpol_orbital_interactions( pose, relevant_residues, struct_id, db_session );
	//report_haro_orbital_interactions( pose, relevant_residues, struct_id, db_session );
	return 0;
}


/// @brief get statistics based upon polar hydrogen to orbital distance/angle
void
OrbitalsFeatures::report_hpol_orbital_interactions(
	Pose const & pose,
	vector1< bool > const &,
	StructureID const struct_id,
	sessionOP db_session
){
	std::string orbita_H_string = "INSERT INTO HPOL_orbital (struct_id, resNum1, orbName1, resNum2, hpolNum2, resName1, orbNum1, resName2, htype2, OrbHdist, cosAOH, cosDHO, chiBAOH, chiBDHO, AOH_angle, DHO_angle, chiBAHD, cosAHD) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);";
	statement orbital_H_statement(basic::database::safely_prepare_statement(orbita_H_string,db_session));
	std::string orbita_orbital_string = "INSERT INTO orbital_orbital (struct_id, resNum1, orbName1, resNum2, orbNum2, resName1, orbNum1, resName2, orbName2, orbOrbdist, cosAOO, cosDOO, chiBAOO, chiBDOO, AOO_angle, DOO_angle, DOA_angle, AOD_angle, chiBAHD, cosAHD) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);";
	statement orbital_orbital_statement(basic::database::safely_prepare_statement(orbita_orbital_string,db_session));
	std::string orbita_Haro_string = "INSERT INTO HARO_orbital (struct_id, resNum1, orbName1, resNum2, haroNum2, resName1, orbNum1, resName2, htype2, orbHdist, cosAOH, cosDHO, chiBAOH, chiBDHO, AOH_angle, DHO_angle, chiBAHD, cosAHD) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);";
	statement orbital_Haro_statement(basic::database::safely_prepare_statement(orbita_Haro_string,db_session));


	for ( Size resNum1 = 1; resNum1 <= pose.size(); ++resNum1 ) {
		Residue  res1 = pose.residue(resNum1);
		Size resNum2(0);
		string orbName1;
		string htype2;
		string resName1=res1.name3();
		string res2name;
		Size orbNum1(0);
		Size hpolNum2(0);
		//Size atm(0);
		Size haroNum2(0);
		Real cosAOH(0);
		Real cosDHO(0);
		Real chiBDHO(0);
		Real chiBAOH(0);
		Real AOH_angle(0);
		Real DHO_angle(0);
		Real chiBAHD(0);
		Real cosAHD(0);
		Size OrbNum2(0);
		string OrbName2;
		//Real OrbOrbdist(0);
		Real cosAOO(0);
		Real cosDOO(0);
		Real chiBAOO(0);
		Real chiBDOO(0);
		Real AOO_angle(0);
		Real DOO_angle(0);
		Real DOA_angle(0);
		Real AOD_angle(0);
		Real OrbHdist(10.1); // min distance used to derive statistics. Should be the shortest distance between an orbital and hydrogen
		bool orb_orb=false;
		bool orb_haro=false;
		bool orb_hpol=false;
		for ( Size res_num2 = 1; res_num2 <= pose.size(); ++res_num2 ) {
			Residue res2 = pose.residue(res_num2);
			if ( resNum1 != res_num2 ) {
				BOOST_FOREACH ( Size const Aindex, res1.atoms_with_orb_index() ) {
					BOOST_FOREACH ( Size const Dindex, res2.atoms_with_orb_index() ) {
						if ( res1.atom_is_backbone(Aindex) || res2.atom_is_backbone(Dindex) ) {
							continue;//just say no to backbone backbone interactions!
						} else {
							BOOST_FOREACH ( Size const Orbindex1, res1.bonded_orbitals(Aindex) ) {
								BOOST_FOREACH ( Size const Orbindex2, res2.bonded_orbitals(Dindex) ) {
									xyzVector<Real> const res1_Orbxyz(res1.orbital_xyz(Orbindex1));
									xyzVector<Real> const res2_Orbxyz(res2.orbital_xyz(Orbindex2));
									Real const container(res1_Orbxyz.distance(res2_Orbxyz));
									//this will actually set OrbHdist to a new low if the container is less than OrbHdist
									//making this type of interaction a pi-pi interaction for cation-pi, not a pi-hpol cation-pi interaction
									if ( container <= OrbHdist ) {
										set_OrbOrb_features_data(res1, res2, Aindex, Dindex, Orbindex1, Orbindex2, res1_Orbxyz, res2_Orbxyz, resNum2,
											orbName1, res2name, OrbName2, orbNum1, OrbNum2, cosAOO, cosDOO, chiBAOO, chiBDOO, AOO_angle,
											DOO_angle, OrbHdist, DOA_angle, AOD_angle, chiBAHD, cosAHD);
										orb_orb=true; //registering with the features reporter that we have a cation_pi interaction
										orb_hpol=false;
										orb_haro=false;
									}
								}
							}
						}
					}
				}
				BOOST_FOREACH ( Size const Aindex, res1.atoms_with_orb_index() ) {
					BOOST_FOREACH ( Size const Orbindex, res1.bonded_orbitals(Aindex) ) {
						xyzVector<Real> const Orbxyz(res1.orbital_xyz(Orbindex));
						BOOST_FOREACH ( Size const Hindex, res2.Haro_index() ) {
							xyzVector<Real> Hxyz(res2.atom(Hindex).xyz());
							Real const container(Orbxyz.distance(Hxyz));
							if ( container <= OrbHdist ) {
								set_OrbH_features_data(res1, res2, Aindex, Hindex, Orbindex, Orbxyz, resNum2, orbName1, htype2, res2name,
									orbNum1, haroNum2, cosAOH, cosDHO, chiBDHO, chiBAOH, AOH_angle, DHO_angle, chiBAHD,cosAHD, OrbHdist );
								orb_haro=true;
								orb_hpol=false;
								orb_orb=false;

							}
						}
					}

				}
				BOOST_FOREACH ( Size const Aindex, res1.atoms_with_orb_index() ) {
					BOOST_FOREACH ( Size const Orbindex, res1.bonded_orbitals(Aindex) ) {
						xyzVector<Real> const Orbxyz(res1.orbital_xyz(Orbindex));

						BOOST_FOREACH ( Size const Hindex, res2.Hpol_index() ) {
							Size Dindex(res2.bonded_neighbor(Hindex)[1]);
							if ( res1.atom_is_backbone(Aindex) && res2.atom_is_backbone(Dindex) ) {
								continue;//do nothing. Dont really want to calculate bb-bb interactions
							} else {
								xyzVector<Real> Hxyz(res2.atom(Hindex).xyz());
								Real const container(Orbxyz.distance(Hxyz));
								if ( container <= OrbHdist ) {
									set_OrbH_features_data(res1, res2, Aindex, Hindex, Orbindex, Orbxyz, resNum2, orbName1, htype2, res2name,
										orbNum1, hpolNum2, cosAOH, cosDHO, chiBDHO, chiBAOH, AOH_angle, DHO_angle, chiBAHD,cosAHD, OrbHdist );
									orb_hpol=true;
									orb_orb=false;
									orb_haro=false;
								}
							}
						}

					}
				}
				//we have to check cation pi interactions. This type of interaction for an arg/aromatic
				//can be between the polar hydrogen of the arg and the pi-orbital of the aromatic ring
				//or it can be between the pi orbital on the arg and the pi-orbital on the aromatic ring
				//this is a special case. Lys does not have a pi-orbital on its nitrogen, therefore, lys
				//can go go on through the if loop that checks distances
				//NOTE this does not take into consideration cation-pi interactions between a ligand
				//something more elegant is needed for that

			}
		}
		if ( OrbHdist <=10.0 && orb_hpol == true ) {
			orbital_H_statement.bind(1,struct_id);
			orbital_H_statement.bind(2, resNum1);
			orbital_H_statement.bind(3, orbName1);
			orbital_H_statement.bind(4, resNum2);
			orbital_H_statement.bind(5, hpolNum2);
			orbital_H_statement.bind(6, resName1);
			orbital_H_statement.bind(7,orbNum1 );
			orbital_H_statement.bind(8, res2name);
			orbital_H_statement.bind(9, htype2);
			orbital_H_statement.bind(10, OrbHdist);
			orbital_H_statement.bind(11, cosAOH);
			orbital_H_statement.bind(12, cosDHO);
			orbital_H_statement.bind(13, chiBAOH);
			orbital_H_statement.bind(14, chiBDHO);
			orbital_H_statement.bind(15, AOH_angle);
			orbital_H_statement.bind(16, DHO_angle);
			orbital_H_statement.bind(17, chiBAHD);
			orbital_H_statement.bind(18, cosAHD);
			basic::database::safely_write_to_database(orbital_H_statement);
		}
		if ( OrbHdist <=10.0 && orb_orb== true ) {
			orbital_orbital_statement.bind(1,struct_id);
			orbital_orbital_statement.bind(2, resNum1);
			orbital_orbital_statement.bind(3, orbName1);
			orbital_orbital_statement.bind(4, resNum2);
			orbital_orbital_statement.bind(5, OrbNum2);
			orbital_orbital_statement.bind(6, resName1);
			orbital_orbital_statement.bind(7, orbNum1);
			orbital_orbital_statement.bind(8, res2name);
			orbital_orbital_statement.bind(9, OrbName2);
			orbital_orbital_statement.bind(10, OrbHdist);
			orbital_orbital_statement.bind(11, cosAOO);
			orbital_orbital_statement.bind(12, cosDOO);
			orbital_orbital_statement.bind(13, chiBAOO);
			orbital_orbital_statement.bind(14, chiBDOO);
			orbital_orbital_statement.bind(15, AOO_angle);
			orbital_orbital_statement.bind(16, DOO_angle);
			orbital_orbital_statement.bind(17, DOA_angle);
			orbital_orbital_statement.bind(18, AOD_angle);
			orbital_orbital_statement.bind(19, chiBAHD);
			orbital_orbital_statement.bind(18, cosAHD);


			basic::database::safely_write_to_database(orbital_orbital_statement);
		}
		if ( OrbHdist <=10.0 && orb_haro == true ) {
			orbital_Haro_statement.bind(1,struct_id);
			orbital_Haro_statement.bind(2, resNum1);
			orbital_Haro_statement.bind(3, orbName1);
			orbital_Haro_statement.bind(4, resNum2);
			orbital_Haro_statement.bind(5, hpolNum2);
			orbital_Haro_statement.bind(6, resName1);
			orbital_Haro_statement.bind(7, orbNum1);
			orbital_Haro_statement.bind(8, res2name);
			orbital_Haro_statement.bind(9, htype2);
			orbital_Haro_statement.bind(10, OrbHdist);
			orbital_Haro_statement.bind(11, cosAOH);
			orbital_Haro_statement.bind(12, cosDHO);
			orbital_Haro_statement.bind(13, chiBAOH);
			orbital_Haro_statement.bind(14, chiBDHO);
			orbital_Haro_statement.bind(15, AOH_angle);
			orbital_Haro_statement.bind(16, DHO_angle);
			orbital_Haro_statement.bind(17, chiBAHD);
			orbital_Haro_statement.bind(18, cosAHD);
			basic::database::safely_write_to_database(orbital_Haro_statement);
		}
	}
}


/// @brief get statistics based upon aromatic hydrogen to orbital distance/angle
void
OrbitalsFeatures::report_haro_orbital_interactions(
	Pose const & /* pose */,
	vector1< bool > const &,
	StructureID const /* struct_id */,
	sessionOP /* db_session */
){
	/* std::string orbita_H_string = "INSERT INTO HARO_orbital (struct_id, resNum1, orbName1, resNum2, haroNum2, resName1, orbNum1, resName2, htype2, orbHdist, cosAOH, cosDHO, chiBAOH, chiBDHO, AOH_angle, DHO_angle, chiBAHD, cosAHD) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);";
	statement orbital_H_statement(basic::database::safely_prepare_statement(orbita_H_string,db_session));


	for(Size resNum1 = 1; resNum1 <= pose.size(); ++resNum1){
	Residue res1 = pose.residue(resNum1);
	Size resNum2(0);
	string orbName1;
	string htype2;
	string res1name;
	string res2name;
	Size orbNum1(0);
	Size haroNum2(0);
	Size atm(0);
	Real cosAOH(0);
	Real cosDHO(0);
	Real chiBDHO(0);
	Real chiBAOH(0);
	Real AOH_angle(0);
	Real DHO_angle(0);
	Real chiBAHD(0);
	Real cosAHD(0);
	Real OrbHdist(10.1); // min distance used to derive statistics. Should be the shortest distance between an orbital and hydrogen
	for(Size res_num2 = 1; res_num2 <= pose.size(); ++res_num2){
	Residue  res2 = pose.residue(res_num2);
	BOOST_FOREACH(Size const Aindex, res1.atoms_with_orb_index()){
	BOOST_FOREACH(Size const Orbindex, res1.bonded_orbitals(Aindex)){
	xyzVector<Real> const Orbxyz(res1.orbital_xyz(Orbindex));

	if(resNum1 != res_num2){
	BOOST_FOREACH(Size const Hindex, res2.Haro_index()){
	xyzVector<Real> Hxyz(res2.atom(Hindex).xyz());
	Real const container(Orbxyz.distance(Hxyz));
	if(container <= OrbHdist){
	set_OrbH_features_data(res1, res2, Aindex, Hindex, Orbindex, Orbxyz, resNum2, orbName1, htype2, res2name,
	orbNum1, haroNum2, cosAOH, cosDHO, chiBDHO, chiBAOH, AOH_angle, DHO_angle, chiBAHD,cosAHD, OrbHdist );

	}
	}
	}
	}

	}
	}
	if(OrbHdist <=10.0){
	orbital_H_statement.bind(1,struct_id);
	orbital_H_statement.bind(2, resNum1);
	orbital_H_statement.bind(3, res1name);
	orbital_H_statement.bind(4,orbNum1 );
	orbital_H_statement.bind(5, orbName1);
	orbital_H_statement.bind(6, resNum2);
	orbital_H_statement.bind(7, res2name);
	orbital_H_statement.bind(8, haroNum2);
	orbital_H_statement.bind(9, htype2);
	orbital_H_statement.bind(10, OrbHdist);
	orbital_H_statement.bind(11, cosAOH);
	orbital_H_statement.bind(12, cosDHO);
	orbital_H_statement.bind(13, chiBAOH);
	orbital_H_statement.bind(14, chiBDHO);
	orbital_H_statement.bind(15, AOH_angle);
	orbital_H_statement.bind(16, DHO_angle);
	orbital_H_statement.bind(17, chiBAHD);
	orbital_H_statement.bind(18, cosAHD);

	basic::database::safely_write_to_database(orbital_H_statement);
	}
	}*/
}


void
OrbitalsFeatures::set_OrbH_features_data(
	Residue & res1,
	Residue & res2,
	Size const Aindex,
	Size const Hindex,
	Size const Orbindex,
	xyzVector<Real> const & Orbxyz,
	Size & resNum2,
	string & orbName1,
	string & htype2,
	string & res2name,
	Size & orbNum1,
	Size & hpolNum2,
	Real & cosAOH,
	Real & cosDHO,
	Real & chiBDHO,
	Real & chiBAOH,
	Real & AOH_angle,
	Real & DHO_angle,
	Real & chiBAHD,
	Real & cosAHD,
	Real & OrbHdist
){
	xyzVector<Real> Hxyz(res2.atom(Hindex).xyz());
	orbName1=res1.orbital_type(Orbindex).name();
	htype2=res2.atom_type(Hindex).name();
	Real const container(Orbxyz.distance(Hxyz));
	xyzVector<Real>  Axyz(res1.atom(Aindex).xyz());
	core::Size Dindex(res2.bonded_neighbor(Hindex)[1]);
	xyzVector<Real> Dxyz(res2.xyz(Dindex));
	AOH_angle = cos_of(Axyz, Orbxyz, Hxyz );
	DHO_angle = cos_of(Dxyz, Hxyz, Orbxyz);
	OrbHdist=container;
	resNum2=res2.seqpos();
	hpolNum2=Hindex;
	std::string res1name=res1.name3();
	orbNum1=Orbindex;
	res2name= res2.name3();

	xyzVector<Real> Bxyz(res2.atom(res2.atom_base(Dindex)).xyz());

	if ( res2name == "PHE" || res2name == "TYR" || res2name == "TRP" || res2name == "HIS" ) {
		res2.update_actcoord();
		Bxyz= res2.actcoord();

	}
	if ( res2name == "ARG" ) {
		Bxyz = res2.atom(res2.atom_index(" CZ ")).xyz();
	}

	xyzVector<Real> AHunit = Orbxyz - Hxyz;
	AHunit.normalize();

	xyzVector<Real> BAunit = Hxyz - Dxyz;
	BAunit.normalize();

	cosDHO = dot( BAunit, AHunit );
	chiBDHO = numeric::dihedral_radians(Bxyz, Dxyz, Hxyz,Orbxyz);

	//BAOH_chi
	Bxyz = res1.atom(res1.atom_base(Aindex)).xyz();

	if ( res1name == "PHE" || res1name == "TYR" || res1name == "TRP" || res1name == "HIS" ) {
		res1.update_actcoord();
		Bxyz= res1.actcoord();
		if ( Orbindex <= 2 ) {
			Bxyz = Axyz;
			Axyz = res1.actcoord();
		}
	}
	if ( res1name == "ARG" ) {
		Bxyz = res1.atom(res1.atom_index(" CZ ")).xyz();
	}
	AHunit = Hxyz - Orbxyz;
	AHunit.normalize();

	BAunit = Orbxyz - Axyz;
	BAunit.normalize();

	cosAOH = dot( BAunit, AHunit );
	chiBAOH = numeric::dihedral_radians(Bxyz,Axyz , Orbxyz,Hxyz);

	AHunit = Dxyz - Hxyz;
	AHunit.normalize();
	BAunit = Hxyz - Axyz;
	BAunit.normalize();

	cosAHD = dot(BAunit, AHunit);
	chiBAHD = numeric::dihedral_radians(Bxyz, Axyz, Hxyz, Dxyz);

}

void
OrbitalsFeatures::set_OrbOrb_features_data(
	Residue & res1,
	Residue & res2,
	Size Aindex,
	Size Dindex,
	Size Orbindex1,
	Size Orbindex2,
	xyzVector<Real> const & Orbxyz1,
	xyzVector<Real> const & Orbxyz2,
	Size & resNum2,
	string & orbName1,
	string & res2name,
	string & OrbName2,
	Size & orbNum1,
	Size & OrbNum2,
	Real & cosAOO,
	Real & cosDOO,
	Real & chiBAOO,
	Real & chiBDOO,
	Real & AOO_angle,
	Real & DOO_angle,
	Real & OrbHdist,
	Real & DOA_angle,
	Real & AOD_angle,
	Real & /* chiBAHD */,
	Real & /* cosAHD */
){
	Real const container(Orbxyz1.distance(Orbxyz2));
	xyzVector<Real> const Axyz(res1.atom(Aindex).xyz());
	xyzVector<Real> Dxyz(res2.xyz(Dindex));
	AOO_angle = cos_of(Axyz, Orbxyz1, Orbxyz2 );
	DOO_angle = cos_of(Dxyz, Orbxyz2, Orbxyz1);
	OrbHdist=container;
	resNum2=res2.seqpos();
	orbName1=res1.orbital_type(Orbindex1).name();
	OrbName2=res2.orbital_type(Orbindex2).name();
	AOD_angle = cos_of(Axyz, Orbxyz1, Dxyz);
	DOA_angle = cos_of(Dxyz, Orbxyz2, Axyz);
	orbNum1=Orbindex1;
	OrbNum2=Orbindex2;
	std::string res1name = res1.name3();
	res2name= res2.name3();

	xyzVector<Real> Bxyz(res2.atom(res2.atom_base(Dindex)).xyz());

	//res2 will never be phe or tyr or trp
	xyzVector<Real> AHunit = Orbxyz1 - Orbxyz2;
	AHunit.normalize();


	xyzVector<Real> BAunit = Orbxyz2 - Dxyz;
	BAunit.normalize();

	cosDOO = dot( BAunit, AHunit );
	chiBDOO = numeric::dihedral_radians(Bxyz, Dxyz, Orbxyz2,Orbxyz1);

	//BAOH_chi
	Bxyz = res1.atom(res1.atom_base(Aindex)).xyz();

	if ( res1name == "PHE" || res1name == "TYR" || res1name == "TRP" || res1name == "HIS" ) {
		res1.update_actcoord();
		Bxyz= res1.actcoord();
	}

	AHunit = Orbxyz2 - Orbxyz1;
	AHunit.normalize();

	BAunit = Orbxyz1 - Axyz;
	BAunit.normalize();

	cosAOO = dot( BAunit, AHunit );
	chiBAOO = numeric::dihedral_radians(Bxyz,Axyz, Orbxyz1,Orbxyz2);
	Bxyz = res1.atom(res1.atom_base(Aindex)).xyz();

	/*
	if(res1name == "PHE" || res1name == "TYR" || res1name == "TRP" || res1name == "HIS"){
	res1.update_actcoord();
	Bxyz= res1.actcoord();
	}

	AHunit = Dxyz - Hxyz;
	AHunit.normalize();
	BAunit = Hxyz - Axyz;
	BAunit.normalize();

	cosAHD = dot(BAunit, AHunit);
	chiBAHD = numeric::dihedral_radians(Bxyz, Axyz, Hxyz, Dxyz);
	*/

}


} // namesapce
} // namespace

