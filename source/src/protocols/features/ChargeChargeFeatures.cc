// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ChargeChargeFeatures.cc
/// @brief  report geometric solvation energy for each polar site to a features database
/// @author Joseph S Harrison

// Unit Headers
#include <protocols/features/ChargeChargeFeatures.hh>

// Platform Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>

//Basic Headers
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>

// Numeric Headers
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/sphericalVector.hh>
#include <core/kinematics/Stub.hh>
#include <numeric/constants.hh>

// Utility Headers
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <string>
#include <sstream>

namespace protocols{
namespace features{

using std::string;
using std::stringstream;
using std::endl;
using core::Size;
using core::Real;
using core::Vector;
using core::Length;
using core::Angle;
using core::Distance;
using core::PointPosition;
using core::kinematics::Stub;
using core::chemical::aa_asp;
using core::chemical::aa_glu;
using core::chemical::aa_lys;
using core::chemical::aa_his;
using core::chemical::aa_arg;
using core::conformation::Residue;
using core::pose::Pose;
using numeric::dihedral_radians;
using numeric::sphericalVector;
using numeric::xyz_to_spherical;
using numeric::constants::r::pi;
using numeric::constants::r::pi_over_2;
using numeric::constants::r::pi_over_180;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;

ChargeChargeFeatures::ChargeChargeFeatures() :
	FeaturesReporter(),
	distance_cutoff_(8)
{}

ChargeChargeFeatures::ChargeChargeFeatures(
	Length distance_cutoff) :
	FeaturesReporter(),
	distance_cutoff_(distance_cutoff)
{}

ChargeChargeFeatures::ChargeChargeFeatures(
	ChargeChargeFeatures const & src ) :
	FeaturesReporter(),
	distance_cutoff_(src.distance_cutoff_)
{}

string
ChargeChargeFeatures::type_name() const { return "ChargeChargeFeatures"; }

void
ChargeChargeFeatures::write_schema_to_db(
	sessionOP db_session
) const {

	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ), false);
	Column q1_site_id("q1_site_id", DbDataTypeOP( new DbInteger() ), false);
	Column q2_site_id("q2_site_id", DbDataTypeOP( new DbInteger() ), false);
	Column q1_charge("q1_charge", DbDataTypeOP( new DbInteger() ), false);
	Column q2_charge("q2_charge", DbDataTypeOP( new DbInteger() ), false);
	Column B1q1q2_angle("B1q1q2_angle", DbDataTypeOP( new DbReal() ), false);
	Column B2q2q1_angle("B2q2q1_angle", DbDataTypeOP( new DbReal() ), false);
	Column q1q2_distance("q1q2_distance", DbDataTypeOP( new DbReal() ), false);
	Column B1q1_torsion("B1q1_torsion", DbDataTypeOP( new DbReal() ), false);
	Column B2q2_torsion("B2q2_torsion", DbDataTypeOP( new DbReal() ), false);

	utility::vector1<Column> pkeys;
	pkeys.push_back(struct_id);
	pkeys.push_back(q1_site_id);
	pkeys.push_back(q2_site_id);


	utility::vector1<Column> fkey_site1_cols;
	fkey_site1_cols.push_back(struct_id);
	fkey_site1_cols.push_back(q1_site_id);

	utility::vector1<std::string> fkey_site1_reference_cols;
	fkey_site1_reference_cols.push_back("struct_id");
	fkey_site1_reference_cols.push_back("site_id");

	utility::vector1<Column> fkey_site2_cols;
	fkey_site2_cols.push_back(struct_id);
	fkey_site2_cols.push_back(q2_site_id);

	utility::vector1<std::string> fkey_site2_reference_cols;
	fkey_site2_reference_cols.push_back("struct_id");
	fkey_site2_reference_cols.push_back("site_id");

	Schema charge_charge_pairs("charge_charge_pairs", PrimaryKey(pkeys));
	charge_charge_pairs.add_column(struct_id);
	charge_charge_pairs.add_column(q1_site_id);
	charge_charge_pairs.add_column(q2_site_id);
	charge_charge_pairs.add_column(q1_charge);
	charge_charge_pairs.add_column(q2_charge);
	charge_charge_pairs.add_column(B1q1q2_angle);
	charge_charge_pairs.add_column(B2q2q1_angle);
	charge_charge_pairs.add_column(q1q2_distance);
	charge_charge_pairs.add_column(B1q1_torsion);
	charge_charge_pairs.add_column(B2q2_torsion);
	charge_charge_pairs.add_foreign_key(
		ForeignKey(fkey_site1_cols, "hbond_sites", fkey_site1_reference_cols, true));
	charge_charge_pairs.add_foreign_key(
		ForeignKey(fkey_site2_cols, "hbond_sites", fkey_site2_reference_cols, true));

	charge_charge_pairs.write(db_session);
}

utility::vector1<std::string>
ChargeChargeFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	dependencies.push_back("HBondFeatures");
	return dependencies;
}

Size
ChargeChargeFeatures::report_features(
	Pose const &,
	vector1< bool > const &,
	StructureID struct_id,
	sessionOP db_session
){

	// locate candidate polar sites from the hbond_sites table

	// Note: these are all pairs of potential salt bridge sites, not
	// just the ones involved in hydrogen bonds
	std::string hbond_string =
		"SELECT\n"
		"	q1.site_id AS q1_id,\n"
		"	q2.site_id AS q2_id,\n"
		"	CASE q1.is_donor WHEN 1 THEN 1 ELSE -1 END AS q1_charge,\n"
		"	CASE q2.is_donor WHEN 1 THEN 1 ELSE -1 END AS q2_charge,\n"
		" q1_xyz.atm_x AS q1_x, q1_xyz.atm_y AS q1_y, q1_xyz.atm_z AS q1_z,\n"
		" q2_xyz.atm_x AS q2_x, q2_xyz.atm_y AS q2_y, q2_xyz.atm_z AS q2_z,\n"
		" q1_xyz.base_x AS B1_x, q1_xyz.base_y AS B1_y, q1_xyz.base_z AS B1_z,\n"
		" q2_xyz.base_x AS B2_x, q2_xyz.base_y AS B2_y, q2_xyz.base_z AS B2_z,\n"
		" q1_xyz.bbase_x AS C1_x, q1_xyz.bbase_y AS C1_y, q1_xyz.bbase_z AS C1_z,\n"
		" q2_xyz.bbase_x AS C2_x, q2_xyz.bbase_y AS C2_y, q2_xyz.bbase_z AS C2_z\n"
		"FROM\n"
		"	hbond_sites AS q1,\n"
		"	hbond_sites AS q2,\n"
		"	hbond_site_atoms AS q1_xyz,\n"
		"	hbond_site_atoms AS q2_xyz\n"
		"WHERE\n"
		"	q1.struct_id = ? AND q2.struct_id = ? AND\n"
		" q1.resNum < q2.resNum AND\n"
		"	(q1.HBChemType = 'hbacc_CXL' OR\n"
		"	q1.HBChemType = 'hbacc_IMD' OR q1.HBChemType = 'hbacc_IME' OR q1.HBChemType = 'hbdon_IMD' OR q1.HBChemType = 'hbdon_IME' OR\n"
		"	q1.HBChemType = 'hbdon_AMO' OR\n"
		"	q1.HBChemType = 'hbdon_GDE' OR q1.HBChemType = 'hbdon_GDH' ) AND\n"
		"	(q2.HBChemType = 'hbacc_CXL' OR\n"
		"	q2.HBChemType = 'hbacc_IMD' OR q2.HBChemType = 'hbacc_IME' OR q2.HBChemType = 'hbdon_IMD' OR q2.HBChemType = 'hbdon_IME' OR\n"
		"	q2.HBChemType = 'hbdon_AMO' OR\n"
		"	q2.HBChemType = 'hbdon_GDE' OR q2.HBChemType = 'hbdon_GDH' ) AND\n"
		"	q1_xyz.struct_id = q1.struct_id AND\n"
		"	q2_xyz.struct_id = q2.struct_id AND\n"
		" q1_xyz.site_id = q1.site_id AND\n"
		" q2_xyz.site_id = q2.site_id;\n";

	statement hbond_stmt(basic::database::safely_prepare_statement(hbond_string,db_session));
	hbond_stmt.bind(1,struct_id);
	hbond_stmt.bind(2,struct_id);

	result res(basic::database::safely_read_from_database(hbond_stmt));
	Size q1_site_id, q2_site_id;
	int q1_charge;
	int q2_charge;
	Length q1_x, q1_y, q1_z;
	Length q2_x, q2_y, q2_z;
	Length B1_x, B1_y, B1_z;
	Length B2_x, B2_y, B2_z;
	Length C1_x, C1_y, C1_z;
	Length C2_x, C2_y, C2_z;

	std::string charge_charge_string = "INSERT INTO charge_charge_pairs (struct_id, q1_site_id, q2_site_id, q1_charge, q2_charge, B1q1q2_angle, B2q2q1_angle, q1q2_distance, B1q1_torsion, B2q2_torsion) VALUES (?,?,?,?,?,?,?,?,?,?);";
	statement charge_charge_statement(basic::database::safely_prepare_statement(charge_charge_string,db_session));


	while(res.next()){
		res >> q1_site_id >> q2_site_id;
		res >> q1_charge >> q2_charge;
		res >> q1_x >> q1_y >> q1_z;
		res >> q2_x >> q2_y >> q2_z;
		res >> B1_x >> B1_y >> B1_z;
		res >> B2_x >> B2_y >> B2_z;
		res >> C1_x >> C1_y >> C1_z;
		res >> C2_x >> C2_y >> C2_z;

		PointPosition const q1(q1_x, q1_y, q1_z);
		PointPosition const q2(q2_x, q2_y, q2_z);
		PointPosition const B1(B1_x, B1_y, B1_z);
		PointPosition const B2(B2_x, B2_y, B2_z);
		PointPosition const C1(C1_x, C1_y, C1_z);
		PointPosition const C2(C2_x, C2_y, C2_z);

		Angle const B1q1q2_angle(angle_of(B1, q1, q2));
		Angle const B2q2q1_angle(angle_of(B2, q2, q1));
		Distance const q1q2_distance = q1.distance(q2);
		Angle const B1q1_torsion(dihedral_radians(C1, B1, q1, q2));
		Angle const B2q2_torsion(dihedral_radians(C2, B2, q2, q1));


		if(q1q2_distance > distance_cutoff_){
			continue;
		}



		charge_charge_statement.bind(1,struct_id);
		charge_charge_statement.bind(2,q1_site_id);
		charge_charge_statement.bind(3,q2_site_id);
		charge_charge_statement.bind(4,q1_charge);
		charge_charge_statement.bind(5,q2_charge);
		charge_charge_statement.bind(6,B1q1q2_angle);
		charge_charge_statement.bind(7,B2q2q1_angle);
		charge_charge_statement.bind(8,q1q2_distance);
		charge_charge_statement.bind(9,B1q1_torsion);
		charge_charge_statement.bind(10,B2q2_torsion);
		basic::database::safely_write_to_database(charge_charge_statement);

	}

	return 0;
}

} //namesapce
} //namespace
