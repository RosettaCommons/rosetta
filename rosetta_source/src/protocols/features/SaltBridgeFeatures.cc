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
using core::chemical::aa_asp;
using core::chemical::aa_glu;
using core::chemical::aa_lys;
using core::chemical::aa_his;
using core::chemical::aa_arg;
using core::conformation::Residue;
using core::pose::Pose;
using numeric::constants::r::pi;
using numeric::dihedral_radians;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;

SaltBridgeFeatures::SaltBridgeFeatures() :
	FeaturesReporter(),
	distance_cutoff_(4)
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
			"	don_site INTEGER,\n"
			"	acc_site INTEGER,\n"
			"	psi REAL,    -- angle around donor group\n"
			"	theta REAL,  -- angle out of donor group plane\n"
			"	rho REAL,    -- distance from center of donor group to acceptor\n"
			"	orbital TEXT,-- syn or anti\n"
			"	FOREIGN KEY (struct_id, don_site)\n"
			"		REFERENCES hbond_sites (struct_id, site_id)\n"
			"		DEFERRABLE INITIALLY DEFERRED,\n"
			"	FOREIGN KEY (struct_id, acc_site)\n"
			"		REFERENCES hbond_sites (struct_id, site_id)\n"
			"		DEFERRABLE INITIALLY DEFERRED,\n"
			"	PRIMARY KEY(struct_id, don_site, acc_site));";
	} else if(db_mode == "mysql") {
		return
			"CREATE TABLE IF NOT EXISTS salt_bridges (\n"
			"	struct_id INTEGER,\n"
			"	don_resNum INTEGER,\n"
			"	acc_resNum INTEGER,\n"
			"	psi REAL,\n"
			"	rho REAL,\n"
			"	theta REAL,\n"
			"	orbital TEXT,\n"
			"	FOREIGN KEY (struct_id, don_site)\n"
			"		REFERENCES hbond_sites (struct_id, site_id)\n"
			"	FOREIGN KEY (struct_id, acc_site)\n"
			"		REFERENCES hbond_sites (struct_id, site_id)\n"
			"	PRIMARY KEY(struct_id, don_site, acc_site));";
	} else {
		return "";
	}
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
	statement stmt = (*db_session) <<
		"SELECT\n"
		"	acc.site_id, acc.resNum, acc.atmNum,\n"
		"	don.site_id, don.resNum, don.atmNum\n"
		"FROM\n"
		"	hbond_sites AS acc, hbond_sites AS don\n"
		"WHERE\n"
		"	don.struct_id = ? AND acc.struct_id = ? AND\n"
		"	(don.HBChemType = 'hbdon_GDH' OR don.HBChemType = 'hbdon_GDE' OR\n"
		"	don.HBChemType = 'hbdon_AMO' OR don.HBChemType = 'hbdon_IMD' OR\n"
		"  don.HBChemType = 'hbdon_IME') AND\n"
		"	(acc.HBChemType = 'hbacc_CXA' OR acc.HBChemType = 'hbacc_CXL');\n"
		<< struct_id << struct_id;

	result res(basic::database::safely_read_from_database(stmt));

	Size acc_site_id, acc_resNum, acc_atmNum;
	Size don_site_id, don_resNum, don_atmNum;
	PointPosition  bb, b, c, o_proj;
	Angle psi, theta;
	while(res.next()){
		res >> acc_site_id >> acc_resNum >> acc_atmNum;
		res >> don_site_id >> don_resNum >> don_atmNum;

		Residue const & d(pose.residue(don_resNum));
		Residue const & a(pose.residue(acc_resNum));
		PointPosition const & n = d.atom(d.atom_base(don_atmNum)).xyz();
		PointPosition const & o = a.atom(acc_atmNum).xyz();
		if(n.distance(o) > distance_cutoff_) continue;


		PointPosition const & ob = a.atom(a.atom_base(acc_atmNum)).xyz();
		PointPosition const & obb = a.atom(a.atom_base(a.atom_base(acc_atmNum))).xyz();
		Angle const hb_chi(dihedral_radians(obb, ob, o, n));
		string const orbital((hb_chi < pi/2 || hb_chi > -pi/2) ? "syn" : "anti");
		switch(d.aa()){
		case aa_his:
			assert(d.atom_name(6) == " CG ");
			b = d.atom(6).xyz();

			assert(d.atom_name(7) == " ND1");
			assert(d.atom_name(10) == " NE2");
			c = (d.atom(7).xyz() + d.atom(10).xyz())/2;

			o_proj = c + dot(b-c, c-o)*(b-c) + dot(n-c, o-c)*(n-c);
			psi = angle_of(b, c, o_proj);
			theta = angle_of(o_proj, c, o);
			break;
		case aa_lys:
			assert(d.atom_name(7) == " CD ");
			bb = d.atom(7).xyz();

			assert(d.atom_name(8) == " CE ");
			b = d.atom(8).xyz();

			assert(d.atom_name(9) == " NZ ");
			c = d.atom(9).xyz();
			psi = dihedral_radians(bb, b, c, o);
			theta = angle_of(b, c, o) - pi/2;
			break;
		case aa_arg:
			assert(d.atom_name(8) == " NE ");
			b = d.atom(8).xyz();

			assert(d.atom_name(9) == " CZ ");
			c = d.atom(9).xyz();
			o_proj = c + dot(b-c, c-o)*(b-c) + dot(n-c, o-c)*(n-c);
			psi = angle_of(b, c, o_proj);
			theta = angle_of(o_proj, c, o);
			break;
		default:
			stringstream err_msg;
			err_msg
				<< "Unrecognized salt bridge donor group: Atom " << d.atom_name(don_atmNum)
				<< "on residue " << don_resNum << "." << endl;
			utility_exit_with_message(err_msg.str());
		}

		Length const rho(c.distance(o));

		statement stmt = (*db_session)
			<< "INSERT INTO salt_bridges VALUES (?,?,?,?,?,?,?)"
			<< struct_id
			<< don_site_id
			<< acc_site_id
			<< psi
			<< theta
			<< rho
			<< orbital;
		basic::database::safely_write_to_database(stmt);
	}

	return 0;
}

} //namesapce
} //namespace
