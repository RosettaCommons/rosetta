// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/SaltBridgeFeatures.cc
/// @brief  report geometric solvation energy for each polar site to a features database
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/SaltBridgeFeatures.hh>

//External
#include <cppdb/frontend.h>

// Platform Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>

// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/database/sql_utils.hh>

// Numeric Headers
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/sphericalVector.hh>
#include <core/kinematics/Stub.hh>
#include <numeric/constants.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>


// Utility Headers
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <string>
#include <sstream>

namespace protocols {
namespace features {

using std::string;
using std::stringstream;
using std::endl;
using core::Size;
using core::Real;
using core::Vector;
using core::Length;
using core::Angle;
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

SaltBridgeFeatures::SaltBridgeFeatures() :
	FeaturesReporter(),
	distance_cutoff_(6)
{}

SaltBridgeFeatures::SaltBridgeFeatures(
	Length distance_cutoff) :
	FeaturesReporter(),
	distance_cutoff_(distance_cutoff)
{}

SaltBridgeFeatures::SaltBridgeFeatures(
	SaltBridgeFeatures const & src ) :
	FeaturesReporter(),
	distance_cutoff_(src.distance_cutoff_)
{}

string
SaltBridgeFeatures::type_name() const { return "SaltBridgeFeatures"; }

void
SaltBridgeFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	write_salt_bridges_table_schema(db_session);
}

void
SaltBridgeFeatures::write_salt_bridges_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column don_resNum("don_resNum", DbDataTypeOP( new DbInteger() ));
	Column acc_id("acc_id", DbDataTypeOP( new DbInteger() ));
	Column psi("psi", DbDataTypeOP( new DbReal() ));
	Column theta("theta", DbDataTypeOP( new DbReal() ));
	Column rho("rho", DbDataTypeOP( new DbReal() ));
	Column orbital("orbital", DbDataTypeOP( new DbText() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(don_resNum);
	primary_key_columns.push_back(acc_id);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns1;
	foreign_key_columns1.push_back(struct_id);
	foreign_key_columns1.push_back(don_resNum);
	vector1< std::string > reference_columns1;
	reference_columns1.push_back("struct_id");
	reference_columns1.push_back("resNum");
	ForeignKey foreign_key1(foreign_key_columns1, "residues", reference_columns1, true);

	Columns foreign_key_columns2;
	foreign_key_columns2.push_back(struct_id);
	foreign_key_columns2.push_back(acc_id);
	vector1< std::string > reference_columns2;
	reference_columns2.push_back("struct_id");
	reference_columns2.push_back("site_id");
	ForeignKey foreign_key2(foreign_key_columns2, "hbond_sites", reference_columns2, true);

	Schema table("salt_bridges", primary_key);
	table.add_foreign_key(foreign_key1);
	table.add_foreign_key(foreign_key2);
	table.add_column(psi);
	table.add_column(theta);
	table.add_column(rho);
	table.add_column(orbital);

	table.write(db_session);
}


utility::vector1<std::string>
SaltBridgeFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	dependencies.push_back("HBondFeatures");
	return dependencies;
}


//Donald JE, Kulp DW, DeGrado WF. Salt bridges: Geometrically
//specific, designable interactions. Proteins: Structure, Function,
//and Bioinformatics. 2010
Size
SaltBridgeFeatures::report_features(
	Pose const & pose,
	vector1< bool > const &,
	StructureID struct_id,
	sessionOP db_session
){

	// locate candidate polar sites from the hbond_sites table

	// Note: these are all pairs of potential salt bridge sites, not
	// just the ones involved in hydrogen bonds
	std::string hbond_string =
		"SELECT\n"
		"\tacc.site_id, acc.resNum, acc.atmNum, don.resNum\n"
		"FROM\n"
		"\thbond_sites AS acc, residues AS don\n"
		"WHERE\n"
		"\tdon.struct_id = ? AND acc.struct_id = ? AND\n"
		" (don.name3 = 'ARG' OR don.name3 = 'LYS' OR don.name3 = 'HIS') AND\n"
		"\t(acc.HBChemType = 'hbacc_CXA' OR acc.HBChemType = 'hbacc_CXL');\n";
	statement hbond_stmt(basic::database::safely_prepare_statement(hbond_string,db_session));
	hbond_stmt.bind(1,struct_id);
	hbond_stmt.bind(2,struct_id);

	result res(basic::database::safely_read_from_database(hbond_stmt));
	Size acc_site_id, acc_resNum, acc_atmNum, don_resNum;
	PointPosition  bb, b, c, n, o_proj;
	Stub don_frame;
	sphericalVector<Real> local_o;
	Angle psi(0), theta(0);
	Length rho(0);
	std::string salt_bridge_string = "INSERT INTO salt_bridges (struct_id, don_resNum, acc_id, psi, theta, rho, orbital) VALUES (?,?,?,?,?,?,?)";
	statement salt_bridge_statement(basic::database::safely_prepare_statement(salt_bridge_string,db_session));


	while ( res.next() ) {
		res >> acc_site_id >> acc_resNum >> acc_atmNum >> don_resNum;

		Residue const & d(pose.residue(don_resNum));
		Residue const & a(pose.residue(acc_resNum));
		PointPosition const & o = a.atom(acc_atmNum).xyz();

		//MJO NOTE: It is necessary to look up the atoms by name rather
		//then using atom numbering directly because residue type patches
		//can insert or delete atoms to modify atom numbering. For
		//instance, the CtermProteinFull patch inserts the OXT atom after
		//the O at position 5 for otherwise unmodified protein backbones.
		switch(d.aa()){
		case aa_lys :
			n = d.atom(d.atom_index(" NZ ")).xyz();
			c = d.atom(d.atom_index(" NZ ")).xyz();
			rho = c.distance(o);
			if ( rho > distance_cutoff_ ) continue;

			bb = d.atom(d.atom_index(" CD ")).xyz();
			b = d.atom(d.atom_index(" CE ")).xyz();

			psi = dihedral_radians(bb, b, c, o);
			theta = angle_of(b, c, o) - pi/2;
			break;
		case aa_his :
			n = d.atom(d.atom_index(" ND1")).xyz();

			c = (n + d.atom(d.atom_index(" NE2")).xyz())/2;
			rho = c.distance(o);
			if ( rho > distance_cutoff_ ) continue;

			b = d.atom(d.atom_index(" CG ")).xyz();
			don_frame.from_four_points(c,c,b,n);
			local_o = xyz_to_spherical(don_frame.global2local(o));
			psi = local_o.phi()*pi_over_180;
			theta = local_o.theta()*pi_over_180 - pi_over_2;
			break;
		case aa_arg :
			c = d.atom(d.atom_index(" CZ ")).xyz();
			rho = c.distance(o);
			if ( rho > distance_cutoff_ ) continue;

			b = d.atom(d.atom_index(" NE ")).xyz();

			{
			PointPosition const & n0(d.atom(d.atom_index(" NE ")).xyz());
			PointPosition const & n1(d.atom(d.atom_index(" NH1")).xyz());
			PointPosition const & n2(d.atom(d.atom_index(" NH2")).xyz());

			n = n0.distance(o) < n1.distance(o) ?
				(n0.distance(o) < n2.distance(o) ? n0 : n2) :
				(n1.distance(o) < n2.distance(o) ? n1 : n2);

			don_frame.from_four_points(c,c,b,n1);
		}

			local_o = xyz_to_spherical(don_frame.global2local(o));
			psi = local_o.phi()*pi_over_180;
			theta = local_o.theta()*pi_over_180 - pi_over_2;
			break;
		default :
			stringstream err_msg;
			err_msg
				<< "Unrecognized salt bridging donor group on rsd: (" << d.name3() << ", " << don_resNum << ")" << endl;
			utility_exit_with_message(err_msg.str());
		}

		PointPosition const & ob = a.atom(a.atom_base(acc_atmNum)).xyz();
		PointPosition const & obb = a.atom(a.atom_base(a.atom_base(acc_atmNum))).xyz();
		Angle const hb_chi(dihedral_radians(obb, ob, o, n));
		string const orbital((hb_chi < pi_over_2 && hb_chi > -pi_over_2) ? "anti" : "syn");

		salt_bridge_statement.bind(1,struct_id);
		salt_bridge_statement.bind(2,don_resNum);
		salt_bridge_statement.bind(3,acc_site_id);
		salt_bridge_statement.bind(4,psi);
		salt_bridge_statement.bind(5,theta);
		salt_bridge_statement.bind(6,rho);
		salt_bridge_statement.bind(7,orbital);
		basic::database::safely_write_to_database(salt_bridge_statement);
	}

	return 0;
}

} //namesapce
} //namespace
