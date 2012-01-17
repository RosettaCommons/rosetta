// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/SaltBridgeFeatures.cc
/// @brief  report geometric solvation energy for each polar site to a features database
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/SaltBridgeFeatures.hh>

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

string
SaltBridgeFeatures::schema() const {
	string db_mode(basic::options::option[basic::options::OptionKeys::inout::database_mode]);

	if(db_mode == "sqlite3")
	{
		return
			"CREATE TABLE IF NOT EXISTS salt_bridges (\n"
			"	struct_id INTEGER,\n"
			"	don_resNum INTEGER,\n"
			"	acc_id INTEGER,\n"
			"	psi REAL,    -- angle around donor group\n"
			"	theta REAL,  -- angle out of donor group plane\n"
			"	rho REAL,    -- distance from center of donor group to acceptor\n"
			"	orbital TEXT,-- syn or anti\n"
			"	FOREIGN KEY (struct_id, don_resNum)\n"
			"		REFERENCES residues (struct_id, resNum)\n"
			"		DEFERRABLE INITIALLY DEFERRED,\n"
			"	FOREIGN KEY (struct_id, acc_id)\n"
			"		REFERENCES hbond_sites (struct_id, site_id)\n"
			"		DEFERRABLE INITIALLY DEFERRED,\n"
			"	PRIMARY KEY(struct_id, don_resNum, acc_id));";
	} else if(db_mode == "mysql") {
		return
			"CREATE TABLE IF NOT EXISTS salt_bridges (\n"
			"	struct_id INTEGER,\n"
			"	don_resNum INTEGER,\n"
			"	acc_id INTEGER,\n"
			"	psi REAL,    -- angle around donor group\n"
			"	theta REAL,  -- angle out of donor group plane\n"
			"	rho REAL,    -- distance from center of donor group to acceptor\n"
			"	orbital TEXT,-- syn or anti\n"
			"	FOREIGN KEY (struct_id, don_site)\n"
			"		REFERENCES hbond_sites (struct_id, site_id)\n"
			"	FOREIGN KEY (struct_id, acc_site)\n"
			"		REFERENCES hbond_sites (struct_id, site_id)\n"
			"	PRIMARY KEY(struct_id, don_site, acc_site));";
	} else {
		return "";
	}
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
	Size struct_id,
	sessionOP db_session
){

	// locate candidate polar sites from the hbond_sites table

	// Note: these are all pairs of potential salt bridge sites, not
	// just the ones involved in hydrogen bonds
	std::string hbond_string =
		"SELECT\n"
		"	acc.site_id, acc.resNum, acc.atmNum, don.resNum\n"
		"FROM\n"
		"	hbond_sites AS acc, residues AS don\n"
		"WHERE\n"
		"	don.struct_id = ? AND acc.struct_id = ? AND\n"
		" (don.name3 = 'ARG' OR don.name3 = 'LYS' OR don.name3 = 'HIS') AND\n"
		"	(acc.HBChemType = 'hbacc_CXA' OR acc.HBChemType = 'hbacc_CXL');\n";
	statement hbond_stmt(basic::database::safely_prepare_statement(hbond_string,db_session));
	hbond_stmt.bind(1,struct_id);
	hbond_stmt.bind(2,struct_id);

	result res(basic::database::safely_read_from_database(hbond_stmt));
	Size acc_site_id, acc_resNum, acc_atmNum, don_resNum;
	PointPosition  bb, b, c, n, o_proj;
	Stub don_frame;
	sphericalVector<Real> local_o;
	Angle psi, theta;
	Length rho;
	std::string salt_bridge_string = "INSERT INTO salt_bridges VALUES (?,?,?,?,?,?,?)";
	statement salt_bridge_statement(basic::database::safely_prepare_statement(salt_bridge_string,db_session));


	while(res.next()){
		res >> acc_site_id >> acc_resNum >> acc_atmNum >> don_resNum;

		Residue const & d(pose.residue(don_resNum));
		Residue const & a(pose.residue(acc_resNum));
		PointPosition const & o = a.atom(acc_atmNum).xyz();

		switch(d.aa()){
		case aa_lys:
			assert(d.atom_name(9) == " NZ ");
			n = d.atom(9).xyz();
			c = d.atom(9).xyz();
			rho = c.distance(o);
			if(rho > distance_cutoff_) continue;

			assert(d.atom_name(7) == " CD ");
			bb = d.atom(7).xyz();

			assert(d.atom_name(8) == " CE ");
			b = d.atom(8).xyz();

			psi = dihedral_radians(bb, b, c, o);
			theta = angle_of(b, c, o) - pi/2;
			break;
		case aa_his:
			assert(d.atom_name(7) == " ND1");
			n = d.atom(7).xyz();

			assert(d.atom_name(10) == " NE2");
			c = (n + d.atom(10).xyz())/2;
			rho = c.distance(o);
			if(rho > distance_cutoff_) continue;

			assert(d.atom_name(6) == " CG ");
			b = d.atom(6).xyz();
			don_frame.from_four_points(c,c,b,n);
			local_o = xyz_to_spherical(don_frame.global2local(o));
			psi = local_o.phi()*pi_over_180;
			theta = local_o.theta()*pi_over_180 - pi_over_2;
			break;
		case aa_arg:
			assert(d.atom_name(9) == " CZ ");
			c = d.atom(9).xyz();
			rho = c.distance(o);
			if(rho > distance_cutoff_) continue;

			{
				assert(d.atom_name(8) == " NE ");
				PointPosition const & n0(d.atom(8).xyz());

				assert(d.atom_name(10) == " NH1");
				PointPosition const & n1(d.atom(10).xyz());

				assert(d.atom_name(11) == " NH2");
				PointPosition const & n2(d.atom(11).xyz());

				n = n0.distance(o) < n1.distance(o) ?
					(n0.distance(o) < n2.distance(o) ? n0 : n2) :
					(n1.distance(o) < n2.distance(o) ? n1 : n2);
			}

			assert(d.atom_name(8) == " NE ");
			b = d.atom(8).xyz();

			don_frame.from_four_points(c,c,b,n);
			local_o = xyz_to_spherical(don_frame.global2local(o));
			psi = local_o.phi()*pi_over_180;
			theta = local_o.theta()*pi_over_180 - pi_over_2;
			break;
		default:
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
