// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ProteinBondGeometryFeatures.cc
/// @brief  report Backbone Torsional Angle features
/// @author Patrick Conway
/// @author Matthew O'Meara

// Unit Headers
#include <core/scoring/methods/CartesianBondedEnergy.hh>
#include <core/scoring/methods/CartBondedParameters.hh>
#include <protocols/features/ProteinBondGeometryFeatures.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>
#include <basic/database/sql_utils.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <numeric/xyz.functions.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/AtomType.hh>
#include <core/kinematics/FoldTree.hh>
#include <basic/basic.hh>
#include <basic/Tracer.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>

#include <core/pose/PDBInfo.hh>

// Platform Headers
#include <core/pose/Pose.hh>

// External Headers
#include <cppdb/frontend.h>
#include <boost/lexical_cast.hpp>

namespace protocols{
namespace features{

using std::string;
using cppdb::statement;
using core::Size;
using core::Real;
using core::conformation::Residue;
using core::pose::Pose;
using utility::vector1;
using utility::sql_database::sessionOP;

static basic::Tracer TR("protocols.features.ProteinBondGeometry");

ProteinBondGeometryFeatures::ProteinBondGeometryFeatures(){
	// if flag _or_ energy method wants a linear potential, make the potential linear - ptc: just the flag here
	linear_bonded_potential_ =
	basic::options::option[ basic::options::OptionKeys::score::linear_bonded_potential ]();

	// initialize databases
	db_ = new core::scoring::methods::IdealParametersDatabase(-1.0,-1.0,-1.0,-1.0,-1.0);
}

ProteinBondGeometryFeatures::ProteinBondGeometryFeatures( ProteinBondGeometryFeatures const & ) :
	FeaturesReporter()
{
}

ProteinBondGeometryFeatures::~ProteinBondGeometryFeatures(){}

string
ProteinBondGeometryFeatures::type_name() const { return "ProteinBondGeometryFeatures"; }

void
ProteinBondGeometryFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	write_bond_intrares_angles_table_schema(db_session);
	write_bond_interres_angles_table_schema(db_session);
	write_bond_intrares_lengths_table_schema(db_session);
	write_bond_interres_lengths_table_schema(db_session);
	write_bond_intrares_torsions_table_schema(db_session);
}

void
ProteinBondGeometryFeatures::write_bond_intrares_angles_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", new DbBigInt(), false);
	Column resNum("resNum", new DbInteger(), false);
	Column cenAtmNum("cenAtmNum", new DbInteger(), false);
	Column outAtm1Num("outAtm1Num", new DbInteger(), false);
	Column outAtm2Num("outAtm2Num", new DbInteger(), false);
	Column cenAtmName("cenAtmName", new DbText(), false);
	Column outAtm1Name("outAtm1Name", new DbText(), false);
	Column outAtm2Name("outAtm2Name", new DbText(), false);
	Column ideal("ideal", new DbReal());
	Column observed("observed", new DbReal());
	Column difference("difference", new DbReal());
	Column energy("energy", new DbReal());

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(resNum);
	primary_key_columns.push_back(cenAtmNum);
	primary_key_columns.push_back(outAtm1Num);
	primary_key_columns.push_back(outAtm2Num);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(struct_id);
	foreign_key_columns.push_back(resNum);
	vector1< std::string > reference_columns;
	reference_columns.push_back("struct_id");
	reference_columns.push_back("resNum");
	ForeignKey foreign_key(foreign_key_columns, "residues", reference_columns, true);

	Schema table("bond_intrares_angles", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(cenAtmName);
	table.add_column(outAtm1Name);
	table.add_column(outAtm2Name);
	table.add_column(ideal);
	table.add_column(observed);
	table.add_column(difference);
	table.add_column(energy);

	table.write(db_session);
}

void
ProteinBondGeometryFeatures::write_bond_interres_angles_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", new DbBigInt(), false);
	Column cenResNum("cenresNum", new DbInteger(), false);
	Column connResNum("connResNum", new DbInteger(), false);
	Column cenAtmNum("cenAtmNum", new DbInteger(), false);
	Column outAtmCenNum("outAtmCenNum", new DbInteger(), false);
	Column outAtmConnNum("outAtmConnNum", new DbInteger(), false);
	Column cenAtmName("cenAtmName", new DbText(), false);
	Column outAtmCenName("outAtmCenName", new DbText(), false);
	Column outAtmConnName("outAtmConnName", new DbText(), false);
	Column ideal("ideal", new DbReal());
	Column observed("observed", new DbReal());
	Column difference("difference", new DbReal());
	Column energy("energy", new DbReal());

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(cenResNum);
	primary_key_columns.push_back(connResNum);
	primary_key_columns.push_back(cenAtmNum);
	primary_key_columns.push_back(outAtmCenNum);
	primary_key_columns.push_back(outAtmConnNum);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(struct_id);
	foreign_key_columns.push_back(cenResNum);
	vector1< std::string > reference_columns;
	reference_columns.push_back("struct_id");
	reference_columns.push_back("resNum");
	ForeignKey foreign_key(foreign_key_columns, "residues", reference_columns, true);

	Schema table("bond_interres_angles", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(cenAtmName);
	table.add_column(outAtmCenName);
	table.add_column(outAtmConnName);
	table.add_column(ideal);
	table.add_column(observed);
	table.add_column(difference);
	table.add_column(energy);

	table.write(db_session);
}

void
ProteinBondGeometryFeatures::write_bond_intrares_lengths_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", new DbBigInt(), false);
	Column resNum("resNum", new DbInteger(), false);
	Column atm1Num("atm1Num", new DbInteger(), false);
	Column atm2Num("atm2Num", new DbInteger(), false);
	Column atm1Name("atm1Name", new DbText(), false);
	Column atm2Name("atm2Name", new DbText(), false);
	Column ideal("ideal", new DbReal());
	Column observed("observed", new DbReal());
	Column difference("difference", new DbReal());
	Column energy("energy", new DbReal());

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(resNum);
	primary_key_columns.push_back(atm1Num);
	primary_key_columns.push_back(atm2Num);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(struct_id);
	foreign_key_columns.push_back(resNum);
	vector1< std::string > reference_columns;
	reference_columns.push_back("struct_id");
	reference_columns.push_back("resNum");
	ForeignKey foreign_key(foreign_key_columns, "residues", reference_columns, true);

	Schema table("bond_intrares_lengths", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(atm1Name);
	table.add_column(atm2Name);
	table.add_column(ideal);
	table.add_column(observed);
	table.add_column(difference);
	table.add_column(energy);

	table.write(db_session);
}

void
ProteinBondGeometryFeatures::write_bond_interres_lengths_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", new DbBigInt(), false);
	Column res1Num("res1Num", new DbInteger(), false);
	Column res2Num("res2Num", new DbInteger(), false);
	Column atm1Num("atm1Num", new DbInteger(), false);
	Column atm2Num("atm2Num", new DbInteger(), false);
	Column atm1Name("atm1Name", new DbText(), false);
	Column atm2Name("atm2Name", new DbText(), false);
	Column ideal("ideal", new DbReal());
	Column observed("observed", new DbReal());
	Column difference("difference", new DbReal());
	Column energy("energy", new DbReal());

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(res1Num);
	primary_key_columns.push_back(res2Num);
	primary_key_columns.push_back(atm1Num);
	primary_key_columns.push_back(atm2Num);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns1;
	foreign_key_columns1.push_back(struct_id);
	foreign_key_columns1.push_back(res1Num);
	vector1< std::string > reference_columns1;
	reference_columns1.push_back("struct_id");
	reference_columns1.push_back("resNum");
	ForeignKey foreign_key1(foreign_key_columns1, "residues", reference_columns1, true);

	Columns foreign_key_columns2;
	foreign_key_columns2.push_back(struct_id);
	foreign_key_columns2.push_back(res2Num);
	vector1< std::string > reference_columns2;
	reference_columns2.push_back("struct_id");
	reference_columns2.push_back("resNum");
	ForeignKey foreign_key2(foreign_key_columns2, "residues", reference_columns2, true);

	Schema table("bond_interres_lengths", primary_key);
	table.add_foreign_key(foreign_key1);
	table.add_foreign_key(foreign_key2);
	table.add_column(atm1Name);
	table.add_column(atm2Name);
	table.add_column(ideal);
	table.add_column(observed);
	table.add_column(difference);
	table.add_column(energy);

	table.write(db_session);
}

void
ProteinBondGeometryFeatures::write_bond_intrares_torsions_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", new DbBigInt(), false);
	Column resNum("resNum", new DbInteger(), false);
	Column atm1Num("atm1Num", new DbInteger(), false);
	Column atm2Num("atm2Num", new DbInteger(), false);
	Column atm3Num("atm3Num", new DbInteger(), false);
	Column atm4Num("atm4Num", new DbInteger(), false);
	Column atm1Name("atm1Name", new DbText(), false);
	Column atm2Name("atm2Name", new DbText(), false);
	Column atm3Name("atm3Name", new DbText(), false);
	Column atm4Name("atm4Name", new DbText(), false);
	Column ideal("ideal", new DbReal(), false);
	Column observed("observed", new DbReal(), false);
	Column difference("difference", new DbReal(), false);
	Column energy("energy", new DbReal(), false);

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(resNum);
	primary_key_columns.push_back(atm1Num);
	primary_key_columns.push_back(atm2Num);
	primary_key_columns.push_back(atm3Num);
	primary_key_columns.push_back(atm4Num);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(struct_id);
	foreign_key_columns.push_back(resNum);
	vector1< std::string > reference_columns;
	reference_columns.push_back("struct_id");
	reference_columns.push_back("resNum");
	ForeignKey foreign_key(foreign_key_columns, "residues", reference_columns, true);

	Schema table("bond_intrares_torsions", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(atm1Name);
	table.add_column(atm2Name);
	table.add_column(atm3Name);
	table.add_column(atm4Name);
	table.add_column(ideal);
	table.add_column(observed);
	table.add_column(difference);
	table.add_column(energy);

	table.write(db_session);
}

utility::vector1<std::string>
ProteinBondGeometryFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	return dependencies;
}

Size
ProteinBondGeometryFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID const struct_id,
	sessionOP db_session
){
	report_intrares_angles( pose, relevant_residues, struct_id, db_session );
	report_interres_angles( pose, relevant_residues, struct_id, db_session );
	report_intrares_lengths( pose, relevant_residues, struct_id, db_session );
	report_interres_lengths( pose, relevant_residues, struct_id, db_session );
	report_intrares_torsions( pose, relevant_residues, struct_id, db_session );
	return 0;
}

void
ProteinBondGeometryFeatures::report_intrares_angles(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID const struct_id,
	sessionOP db_session
){
	std::string statement_string ="INSERT INTO bond_intrares_angles (struct_id, resNum, cenAtmNum, outAtm1Num, outAtm2Num, cenAtmName, outAtm1Name, outAtm2Name, ideal, observed, difference, energy) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));

	Real energy_angle = 0;

	for (Size i = 1; i <= pose.total_residue(); ++i) {
		if(!check_relevant_residues(relevant_residues, i)) continue;

		Residue const & rsd = pose.residue(i);
		if(!rsd.is_protein()) continue;

		//following code ripped off from core/scoring/methods/CartesianBondedEnergy.cc
		//energy computations are not up to date with current cart_bonded - the rest is ok

		// get residue type
		core::chemical::ResidueType const & rsd_type = rsd.type();

		// for each angle in the residue
		for ( Size bondang = 1; bondang <= rsd_type.num_bondangles(); ++bondang ) {
			// get ResidueType ints
			Size rt1 = ( rsd_type.bondangle( bondang ) ).key1();
			Size rt2 = ( rsd_type.bondangle( bondang ) ).key2();
			Size rt3 = ( rsd_type.bondangle( bondang ) ).key3();

			// check for vrt
			//if ( rsd_type.atom_type(rt1).is_virtual()
			//       || rsd_type.atom_type(rt2).is_virtual()
			//       || rsd_type.atom_type(rt3).is_virtual() )
			if ( rsd_type.aa() == core::chemical::aa_vrt)
				continue;

			// lookup Ktheta and theta0
			Real Ktheta, theta0;
			db_->lookup_angle_legacy( pose, rsd, rt1, rt2, rt3, Ktheta, theta0 );
			if (Ktheta == 0.0) continue;

			// get angle
			Real const angle = numeric::angle_radians(
														rsd.atom( rt1 ).xyz(),
														rsd.atom( rt2 ).xyz(),
														rsd.atom( rt3 ).xyz() );

			if (linear_bonded_potential_ && std::fabs(angle - theta0)>1) {
				energy_angle = 0.5*Ktheta*std::fabs(angle-theta0);
				//TR << "intrares_angles - linear_bonded energy: " << energy_angle << std::endl;
			} else {
				energy_angle = 0.5*Ktheta*(angle-theta0) * (angle-theta0);
				//TR << "intrares_angles - energy: " << energy_angle << std::endl;
				/*TR.Debug << pose.pdb_info()->name() << " seqpos: " << rsd.seqpos() << " pdbpos: " << pose.pdb_info()->number(rsd.seqpos()) << " intrares angle: " <<
 				rsd_type.name() << " : " <<
				rsd.atom_name( rt1 ) << " , " << rsd.atom_name( rt2 ) << " , " <<
				rsd.atom_name( rt3 ) << "   " << angle << "  " << theta0 << "     " <<
				Ktheta << " " << 0.5*Ktheta*(angle-theta0) * (angle-theta0) << std::endl;*/
			}

			const std::string tmp = boost::lexical_cast<std::string>(struct_id);

			//report results here
			stmt.bind(1,struct_id);
			stmt.bind(2,i);
			stmt.bind(3,rt2);
			stmt.bind(4,rt1);
			stmt.bind(5,rt3);
			stmt.bind(6,rsd.atom_name( rt2 ));
			stmt.bind(7,rsd.atom_name( rt1 ));
			stmt.bind(8,rsd.atom_name( rt3 ));
			stmt.bind(9,theta0);
			stmt.bind(10,angle);
			stmt.bind(11,angle-theta0);
			stmt.bind(12,energy_angle);
			basic::database::safely_write_to_database(stmt);
		}
	}
}

void
ProteinBondGeometryFeatures::report_interres_angles(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID const struct_id,
	sessionOP db_session
){
	std::string statement_string ="INSERT INTO bond_interres_angles (struct_id, cenresNum, connResNum, cenAtmNum, outAtmCenNum, outAtmConnNum, cenAtmName, outAtmCenName, outAtmConnName, ideal, observed, difference, energy) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));

	for (Size i = 1; i <= pose.total_residue(); ++i) {
		Residue const & rsd1 = pose.residue(i);
		if(!rsd1.is_protein()) continue;

		for (Size j = i+1; j <= pose.total_residue(); ++j) {
			if(!check_relevant_residues(relevant_residues, i, j)) continue;
			Residue const & rsd2 = pose.residue(j);
			if(!rsd2.is_protein()) continue;

			//following code ripped off from core/scoring/methods/CartesianBondedEnergy.cc
		  //energy computations are not up to date with current cart_bonded - the rest is ok

			// bail out if the residues aren't bonded
			if (!rsd1.is_bonded(rsd2)) continue;

			//fpd chainbreak variants also mess things up
			//fpd check for chainbreaks
			if ( pose.fold_tree().is_cutpoint( std::min( rsd1.seqpos(), rsd2.seqpos() ) ) ) continue;

			// get residue types
			core::chemical::ResidueType const & rsd1_type = rsd1.type();
			core::chemical::ResidueType const & rsd2_type = rsd2.type();

			utility::vector1< Size > const & r1_resconn_ids( rsd1.connections_to_residue( rsd2 ) );

			for ( Size ii = 1; ii <= r1_resconn_ids.size(); ++ii ) {

				Size const resconn_id1( r1_resconn_ids[ii] );
				Size const resconn_id2( rsd1.residue_connection_conn_id( resconn_id1 ) );

				Size const resconn_atomno1( rsd1.residue_connection( resconn_id1 ).atomno() );
				Size const resconn_atomno2( rsd2.residue_connection( resconn_id2 ).atomno() );

				/// compute the bond-angle energies from pairs of atoms within-1 bond on rsd1 with
				/// the the connection atom on rsd2.
				utility::vector1< core::chemical::two_atom_set > const & rsd1_atoms_wi1_bond_of_ii(
																							 rsd1_type.atoms_within_one_bond_of_a_residue_connection( resconn_id1 ));
				for ( Size jj = 1; jj <= rsd1_atoms_wi1_bond_of_ii.size(); ++jj ) {
					assert( rsd1_atoms_wi1_bond_of_ii[ jj ].key1() == resconn_atomno1 );
					Size const res1_lower_atomno = rsd1_atoms_wi1_bond_of_ii[ jj ].key2();

					Real const angle = numeric::angle_radians(
																rsd1.atom( res1_lower_atomno ).xyz(),
																rsd1.atom( resconn_atomno1 ).xyz(),
																rsd2.atom( resconn_atomno2 ).xyz() );

					// lookup Ktheta and theta0
					Real Ktheta, theta0;
					db_->lookup_angle_legacy( pose, rsd1, res1_lower_atomno, resconn_atomno1, -resconn_id1, Ktheta, theta0 );

					if (Ktheta == 0.0) continue;

					// accumulate the energy
					Real energy_angle = 0;		//ptc - don't accumulate, report each angle on it's own
					if (linear_bonded_potential_ && std::fabs(angle-theta0)>1) {
						energy_angle += 0.5*Ktheta*std::fabs(angle-theta0);
					} else {
						energy_angle += 0.5*Ktheta*(angle-theta0) * (angle-theta0);
					}

					//report results here
					stmt.bind(1,struct_id);
					stmt.bind(2,i);
					stmt.bind(3,j);
					stmt.bind(4,resconn_atomno1);
					stmt.bind(5,res1_lower_atomno);
					stmt.bind(6, resconn_atomno2);
					stmt.bind(7, rsd1.atom_name( resconn_atomno1 ));
					stmt.bind(8, rsd1.atom_name( res1_lower_atomno ));
					stmt.bind(9, rsd2.atom_name( resconn_atomno2 ));
					stmt.bind(10,theta0);
					stmt.bind(11,angle);
					stmt.bind(12,angle-theta0);
					stmt.bind(13,energy_angle);
					basic::database::safely_write_to_database(stmt);
				}

				/// compute the bond-angle energies from pairs of atoms within-1 bond on rsd2 with
				/// the the connection atom on rsd1.
				utility::vector1< core::chemical::two_atom_set > const & rsd2_atoms_wi1_bond_of_ii(
																							 rsd2_type.atoms_within_one_bond_of_a_residue_connection( resconn_id2 ));
				for ( Size jj = 1; jj <= rsd2_atoms_wi1_bond_of_ii.size(); ++jj ) {
					assert( rsd2_atoms_wi1_bond_of_ii[ jj ].key1() == resconn_atomno2 );
					Size const res2_lower_atomno = rsd2_atoms_wi1_bond_of_ii[ jj ].key2();

					// lookup Ktheta and theta0
					Real Ktheta, theta0;
					db_->lookup_angle_legacy( pose, rsd2, res2_lower_atomno, resconn_atomno2, -resconn_id2, Ktheta, theta0 );

					if (Ktheta == 0.0) continue;
					Real const angle = numeric::angle_radians(
																rsd2.atom( res2_lower_atomno ).xyz(),
																rsd2.atom( resconn_atomno2 ).xyz(),
																rsd1.atom( resconn_atomno1 ).xyz() );

					// accumulate the energy
					Real energy_angle = 0;		//ptc - don't accumulate, report each angle on it's own
					if (linear_bonded_potential_ && std::fabs(angle-theta0)>1) {
						energy_angle += 0.5*Ktheta*std::fabs(angle-theta0);
					} else {
						energy_angle += 0.5*Ktheta*(angle-theta0) * (angle-theta0);
					}

					//report results here
					stmt.bind(1,struct_id);
					stmt.bind(2,j);
					stmt.bind(3,i);
					stmt.bind(4,resconn_atomno2);
					stmt.bind(5,res2_lower_atomno);
					stmt.bind(6, resconn_atomno1);
					stmt.bind(7, rsd2.atom_name( resconn_atomno2 ));
					stmt.bind(8, rsd2.atom_name( res2_lower_atomno ));
					stmt.bind(9, rsd1.atom_name( resconn_atomno1 ));
					stmt.bind(10,theta0);
					stmt.bind(11,angle);
					stmt.bind(12,angle-theta0);
					stmt.bind(13,energy_angle);
					basic::database::safely_write_to_database(stmt);
				}
			}
		}
	}
}

void
ProteinBondGeometryFeatures::report_intrares_lengths(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID const struct_id,
	sessionOP db_session
){
	std::string statement_string ="INSERT INTO bond_intrares_lengths (struct_id, resNum, atm1Num, atm2Num, atm1Name, atm2Name, ideal, observed, difference, energy) VALUES (?,?,?,?,?,?,?,?,?,?)";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));

	for (Size i = 1; i <= pose.total_residue(); ++i) {
		if(!check_relevant_residues(relevant_residues, i)) continue;

		Residue const & rsd = pose.residue(i);
		if(!rsd.is_protein()) continue;

		//following code ripped off from core/scoring/methods/CartesianBondedEnergy.cc
		//energy computations are not up to date with current cart_bonded - the rest is ok

		core::chemical::ResidueType const & rsd_type = rsd.type();

		// for each bond in the residue
		// for each bonded atom
		for (Size atm_i=1; atm_i<=rsd_type.natoms(); ++atm_i) {
			core::chemical::AtomIndices atm_nbrs = rsd_type.nbrs( atm_i );
			for (Size j=1; j<=atm_nbrs.size(); ++j) {
				Size atm_j = atm_nbrs[j];
				if ( atm_i<atm_j ) { // only score each bond once -- use restype index to define ordering
					// check for vrt
					//if ( rsd_type.atom_type(atm_i).is_virtual() || rsd_type.atom_type(atm_j).is_virtual() )
					if ( rsd_type.aa() == core::chemical::aa_vrt)
						continue;

					// lookup Ktheta and theta0
					Real Kd, d0;
					db_->lookup_length_legacy( pose, rsd, atm_i, atm_j, Kd, d0 );
					if (Kd == 0.0) continue;

					Real const d = ( rsd.atom( atm_i ).xyz()-rsd.atom( atm_j ).xyz() ).length();

					// accumulate the energy
					Real energy_length = 0;		//ptc - don't accumulate, report each length on it's own
					if (linear_bonded_potential_ && std::fabs(d - d0)>1) {
						energy_length += 0.5*Kd*std::fabs(d-d0);
					} else {
						energy_length += 0.5*Kd*(d-d0)*(d-d0);
					}

					//report results here
					stmt.bind(1,struct_id);
					stmt.bind(2,i);
					stmt.bind(3,atm_i);
					stmt.bind(4,atm_j);
					stmt.bind(5, rsd.atom_name( atm_i ));
					stmt.bind(6, rsd.atom_name( atm_j ));
					stmt.bind(7,d0);
					stmt.bind(8,d);
					stmt.bind(9,d-d0);
					stmt.bind(10,energy_length);
					basic::database::safely_write_to_database(stmt);
				}
			}
		}
	}
}

void
ProteinBondGeometryFeatures::report_interres_lengths(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID const struct_id,
	sessionOP db_session
){
	std::string statement_string ="INSERT INTO bond_interres_lengths (struct_id, res1Num, res2Num, atm1Num, atm2Num, atm1Name, atm2Name, ideal, observed, difference, energy) VALUES (?,?,?,?,?,?,?,?,?,?,?)";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));


	for (Size i = 1; i <= pose.total_residue(); ++i) {
		Residue const & rsd1 = pose.residue(i);
		if(!rsd1.is_protein()) continue;

		for (Size j = i+1; j <= pose.total_residue(); ++j) {
			if(!check_relevant_residues(relevant_residues, i, j)) continue;
			Residue const & rsd2 = pose.residue(j);
			if(!rsd2.is_protein()) continue;

			//following code ripped off from core/scoring/methods/CartesianBondedEnergy.cc
		  //energy computations are not up to date with current cart_bonded - the rest is ok

			// bail out if the residues aren't bonded
			if (!rsd1.is_bonded(rsd2)) continue;

			//fpd chainbreak variants also mess things up
			//fpd check for chainbreaks
			if ( pose.fold_tree().is_cutpoint( std::min( rsd1.seqpos(), rsd2.seqpos() ) ) ) continue;

			utility::vector1< Size > const & r1_resconn_ids( rsd1.connections_to_residue( rsd2 ) );

			for ( Size ii = 1; ii <= r1_resconn_ids.size(); ++ii ) {

				Size const resconn_id1( r1_resconn_ids[ii] );
				Size const resconn_id2( rsd1.residue_connection_conn_id( resconn_id1 ) );

				Size const resconn_atomno1( rsd1.residue_connection( resconn_id1 ).atomno() );
				Size const resconn_atomno2( rsd2.residue_connection( resconn_id2 ).atomno() );


				/// finally, compute the bondlength across the interface
				Real length =
				( rsd2.atom( resconn_atomno2 ).xyz() - rsd1.atom( resconn_atomno1 ).xyz() ).length();

				// lookup Ktheta and theta0
				Real Kd, d0;
				db_->lookup_length_legacy( pose, rsd1, resconn_atomno1, -resconn_id1, Kd, d0 );

				// accumulate the energy
				Real energy_length = 0;			//ptc - dont accumulate energy, report each length on it's own.
				if (linear_bonded_potential_ && std::fabs(length-d0)>1) {
					energy_length += 0.5*Kd*std::fabs(length-d0);
				} else {
					energy_length += 0.5*Kd*(length-d0)*(length-d0);
				}

				//report results here
				stmt.bind(1,struct_id);
				stmt.bind(2,i);
				stmt.bind(3,j);
				stmt.bind(4,resconn_atomno1);
				stmt.bind(5,resconn_atomno2);
				stmt.bind(6, rsd1.atom_name( resconn_atomno1 ));
				stmt.bind(7, rsd2.atom_name( resconn_atomno2 ));
				stmt.bind(8,d0);
				stmt.bind(9,length);
				stmt.bind(10,length-d0);
				stmt.bind(11,energy_length);
				basic::database::safely_write_to_database(stmt);
			}
		}
	}
}

void
ProteinBondGeometryFeatures::report_intrares_torsions(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID const struct_id,
	sessionOP db_session
){
	std::string statement_string ="INSERT INTO bond_intrares_torsions (struct_id, resNum, atm1Num, atm2Num, atm3Num, atm4Num, atm1Name, atm2Name, atm3Name, atm4Name, ideal, observed, difference, energy) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));

	for (Size i = 1; i <= pose.total_residue(); ++i) {
		if(!check_relevant_residues(relevant_residues, i)) continue;

		Residue const & rsd = pose.residue(i);
		if(!rsd.is_protein()) continue;

		//following code ripped off from core/scoring/methods/CartesianBondedEnergy.cc
		//energy computations are not up to date with current cart_bonded - the rest is ok

		core::chemical::ResidueType const & rsd_type = rsd.type();

		// for each torsion _that doesn't correspond to a DOF_ID in the pose_
		for ( Size dihe = 1; dihe <= rsd_type.ndihe(); ++dihe ){
			// get ResidueType ints
			int rt1 = ( rsd_type.dihedral( dihe ) ).key1();
			int rt2 = ( rsd_type.dihedral( dihe ) ).key2();
			int rt3 = ( rsd_type.dihedral( dihe ) ).key3();
			int rt4 = ( rsd_type.dihedral( dihe ) ).key4();

			// lookup Ktheta and theta0
			Real Kphi, phi0, phi_step;
			db_->lookup_torsion_legacy( rsd.type(), rt1, rt2, rt3, rt4, Kphi, phi0, phi_step );
			if (Kphi == 0.0) continue;

			// get angle
			Real angle = numeric::dihedral_radians
			( rsd.atom( rt1 ).xyz(), rsd.atom( rt2 ).xyz(),
			 rsd.atom( rt3 ).xyz(), rsd.atom( rt4 ).xyz() );

			// accumulate the energy
			Real energy_torsion = 0;			//ptc - dont accumulate energy, report each torsion on it's own
			Real del_phi = basic::subtract_radian_angles(angle, phi0);
			if (phi_step>0) del_phi = basic::periodic_range( del_phi, phi_step );

			if (linear_bonded_potential_ && std::fabs(del_phi)>1)
				energy_torsion += 0.5*Kphi*std::fabs(del_phi);
			else
				energy_torsion += 0.5*Kphi*del_phi*del_phi;

			//report results here
			stmt.bind(1,struct_id);
			stmt.bind(2,i);
			stmt.bind(3,rt1);
			stmt.bind(4,rt2);
			stmt.bind(5,rt3);
			stmt.bind(6,rt4);
			stmt.bind(7,rsd.atom_name( rt1 ) );
			stmt.bind(8,rsd.atom_name( rt2 ) );
			stmt.bind(9,rsd.atom_name( rt3 ) );
			stmt.bind(10,rsd.atom_name( rt4 ) );
			stmt.bind(11,phi0);
			stmt.bind(12,angle);
			stmt.bind(13,del_phi);
			stmt.bind(14,energy_torsion);
			basic::database::safely_write_to_database(stmt);
		}
	}
}

} // namesapce
} // namespace
