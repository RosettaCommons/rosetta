// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ResidueConformationFeatures.cc
/// @brief  report idealized torsional DOFs Statistics Scientific Benchmark
/// @author Matthew O'Meara
/// @author Sam DeLuca

// Unit Headers
#include <protocols/features/ResidueConformationFeatures.hh>

//External

// Project Headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/types.hh>
#include <core/id/AtomID.hh>

// Utility Headers
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/tools/make_vector.hh>
#include <protocols/features/util.hh>

//Basic Headers
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/Constraint.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/database/insert_statement_generator/InsertGenerator.hh>
#include <basic/database/insert_statement_generator/RowData.hh>

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
using basic::database::table_exists;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;
using basic::database::insert_statement_generator::InsertGenerator;
using basic::database::insert_statement_generator::RowDataBaseOP;
using basic::database::insert_statement_generator::RowData;

ResidueConformationFeatures::ResidueConformationFeatures()
{
	compact_residue_schema_ = basic::options::option[basic::options::OptionKeys::inout::dbms::use_compact_residue_schema]();
}

ResidueConformationFeatures::ResidueConformationFeatures(
		ResidueConformationFeatures const & ) : FeaturesReporter()
{
	compact_residue_schema_ = basic::options::option[basic::options::OptionKeys::inout::dbms::use_compact_residue_schema]();
}

string
ResidueConformationFeatures::type_name() const {
	return "ResidueConformationFeatures";
}

void
ResidueConformationFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const{
	using namespace basic::database::schema_generator;

	//******nonprotein_residue_conformation******//
	Column struct_id("struct_id", new DbBigInt(), false);
	Column seqpos("seqpos", new DbInteger(), false);
	Column phi("phi", new DbDouble(), false);
	Column psi("psi", new DbDouble(), false);
	Column omega("omega", new DbDouble(), false);

	utility::vector1<Column> non_prot_res_pkeys;
	non_prot_res_pkeys.push_back(struct_id);
	non_prot_res_pkeys.push_back(seqpos);



	utility::vector1<Column> fkey_cols;
	fkey_cols.push_back(struct_id);
	fkey_cols.push_back(seqpos);

	utility::vector1<std::string> fkey_reference_cols;
	fkey_reference_cols.push_back("struct_id");
	fkey_reference_cols.push_back("resNum");

	Schema nonprotein_residue_conformation("nonprotein_residue_conformation", PrimaryKey(non_prot_res_pkeys));
	nonprotein_residue_conformation.add_column(struct_id);
	nonprotein_residue_conformation.add_column(seqpos);
	nonprotein_residue_conformation.add_column(phi);
	nonprotein_residue_conformation.add_column(psi);
	nonprotein_residue_conformation.add_column(omega);
	nonprotein_residue_conformation.add_foreign_key(ForeignKey(fkey_cols, "residues", fkey_reference_cols, true));

	nonprotein_residue_conformation.write(db_session);

	//******nonprotein_residue_angles******//
	Column chinum("chinum", new DbInteger(), false);
	Column chiangle("chiangle", new DbDouble(), false);

	utility::vector1<Column> non_prot_res_angle_keys;
	non_prot_res_angle_keys.push_back(struct_id);
	non_prot_res_angle_keys.push_back(seqpos);
	non_prot_res_angle_keys.push_back(chinum);

	Schema nonprotein_residue_angles("nonprotein_residue_angles", PrimaryKey(non_prot_res_angle_keys));
	nonprotein_residue_angles.add_column(struct_id);
	nonprotein_residue_angles.add_column(seqpos);
	nonprotein_residue_angles.add_column(chinum);
	nonprotein_residue_angles.add_column(chiangle);
	nonprotein_residue_angles.add_foreign_key(ForeignKey(fkey_cols, "residues", fkey_reference_cols, true));

	nonprotein_residue_angles.write(db_session);

	if(compact_residue_schema_)
	{
		//******compact_residue_atom_coords*****//
		Column coord_data("coord_data",new DbText());
		Column atom_count("atom_count",new DbInteger(),false);
		utility::vector1<Column> compact_res_atom_coords_pkeys;
		compact_res_atom_coords_pkeys.push_back(struct_id);
		compact_res_atom_coords_pkeys.push_back(seqpos);

		Schema compact_residue_atom_coords("compact_residue_atom_coords",PrimaryKey(compact_res_atom_coords_pkeys));
		compact_residue_atom_coords.add_column(struct_id);
		compact_residue_atom_coords.add_column(seqpos);
		compact_residue_atom_coords.add_column(atom_count);
		compact_residue_atom_coords.add_column(coord_data);
		compact_residue_atom_coords.add_foreign_key(ForeignKey(fkey_cols, "residues", fkey_reference_cols, true));
		compact_residue_atom_coords.write(db_session);

	}else
	{
		//******residue_atom_coords******//
		Column atomno("atomno", new DbInteger(), false);
		Column x("x", new DbDouble(), false);
		Column y("y", new DbDouble(), false);
		Column z("z", new DbDouble(), false);

		utility::vector1<Column> res_atm_coords_pkeys;
		res_atm_coords_pkeys.push_back(struct_id);
		res_atm_coords_pkeys.push_back(seqpos);
		res_atm_coords_pkeys.push_back(atomno);

		Schema residue_atom_coords("residue_atom_coords", PrimaryKey(res_atm_coords_pkeys));
		residue_atom_coords.add_column(struct_id);
		residue_atom_coords.add_column(seqpos);
		residue_atom_coords.add_column(atomno);
		residue_atom_coords.add_column(x);
		residue_atom_coords.add_column(y);
		residue_atom_coords.add_column(z);
		residue_atom_coords.add_foreign_key(ForeignKey(fkey_cols, "residues", fkey_reference_cols, true));

		residue_atom_coords.write(db_session);
	}
}

utility::vector1<std::string>
ResidueConformationFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	return dependencies;
}

Size
ResidueConformationFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID struct_id,
	sessionOP db_session
){
	//check to see if this structure is ideal

	bool ideal = true;
	if(!basic::options::option[basic::options::OptionKeys::out::file::force_nonideal_structure]())
	{
		core::conformation::Conformation const & conformation(pose.conformation());
		for(core::Size resn=1; resn <= pose.n_residue();++resn){
			if(!check_relevant_residues(relevant_residues, resn)) continue;
			bool residue_status(core::conformation::is_ideal_position(resn,conformation));
			if(!residue_status){
				ideal = false;
				break;
			}
		}
	}else
	{
		ideal = false;
	}

	InsertGenerator conformation_insert("nonprotein_residue_conformation");
	conformation_insert.add_column("struct_id");
	conformation_insert.add_column("seqpos");
	conformation_insert.add_column("phi");
	conformation_insert.add_column("psi");
	conformation_insert.add_column("omega");

	InsertGenerator angle_insert("nonprotein_residue_angles");
	angle_insert.add_column("struct_id");
	angle_insert.add_column("seqpos");
	angle_insert.add_column("chinum");
	angle_insert.add_column("chiangle");

	//We only need one of these but the scoping gets easier to deal with if we just make both
	//InsertGenerators are really lightweight and everything else about database IO is slower than this function anyways
	//This entire function is in desperate need of some method extraction, if you're feeling like doing a good deed, this would be one
	InsertGenerator atom_insert("residue_atom_coords");
	atom_insert.add_column("struct_id");
	atom_insert.add_column("seqpos");
	atom_insert.add_column("atomno");
	atom_insert.add_column("x");
	atom_insert.add_column("y");
	atom_insert.add_column("z");

	InsertGenerator compact_residue_insert("compact_residue_atom_coords");
	compact_residue_insert.add_column("struct_id");
	compact_residue_insert.add_column("seqpos");
	compact_residue_insert.add_column("atom_count");
	compact_residue_insert.add_column("coord_data");

	RowDataBaseOP struct_id_data = new RowData<StructureID>("struct_id",struct_id);
	for (Size i = 1; i <= pose.total_residue(); ++i) {
		if(!check_relevant_residues(relevant_residues, i)) continue;

		Residue const & resi = pose.residue(i);
		if(resi.aa() <= num_canonical_aas){
			continue;
		}
		//runtime_assert(resi.aa() <= num_canonical_aas);
		Real phi  (0.0);
		Real psi  (0.0);
		Real omega(0.0);
		if(!resi.is_ligand()){
			phi  = resi.mainchain_torsion(1);
			psi  = resi.mainchain_torsion(2);
			omega= resi.mainchain_torsion(3);
		}

		RowDataBaseOP seqpos_data = new RowData<Size>("seqpos",i);
		RowDataBaseOP phi_data = new RowData<Real>("phi",phi);
		RowDataBaseOP psi_data = new RowData<Real>("psi",psi);
		RowDataBaseOP omega_data = new RowData<Real>("omega",omega);

		conformation_insert.add_row(
			utility::tools::make_vector(struct_id_data,seqpos_data,phi_data,psi_data,omega_data));



		for(core::Size chi_num = 1; chi_num <= resi.nchi();++chi_num){
			core::Real chi_angle = resi.chi(chi_num);

			RowDataBaseOP chinum_data = new RowData<Size>("chinum",chi_num);
			RowDataBaseOP chiangle_data = new RowData<Real>("chiangle",chi_angle);
			angle_insert.add_row(utility::tools::make_vector(
				struct_id_data,seqpos_data,chinum_data,chiangle_data));

		}
		if(!ideal || resi.is_ligand()){ // always store coords for a ligand

			if(compact_residue_schema_)
			{

				std::string residue_data_string(serialize_residue_xyz_coords(resi));
				RowDataBaseOP atom_count = new RowData<core::Size>("atom_count",resi.natoms());
				RowDataBaseOP residue_data = new RowData<std::string>("coord_data",residue_data_string);
				compact_residue_insert.add_row(utility::tools::make_vector(struct_id_data,atom_count,seqpos_data,residue_data));
			}else
			{
				for(Size atom = 1; atom <= resi.natoms(); ++atom){
					core::Vector coords = resi.xyz(atom);

					RowDataBaseOP atom_data = new RowData<Size>("atomno",atom);
					RowDataBaseOP x_data = new RowData<Real>("x",coords.x());
					RowDataBaseOP y_data = new RowData<Real>("y",coords.y());
					RowDataBaseOP z_data = new RowData<Real>("z",coords.z());

					atom_insert.add_row(utility::tools::make_vector(
						struct_id_data,seqpos_data,atom_data,x_data,y_data,z_data));
				}
			}
		}
	}

	conformation_insert.write_to_database(db_session);
	angle_insert.write_to_database(db_session);
	if(compact_residue_schema_)
	{
		compact_residue_insert.write_to_database(db_session);
	}else
	{
		atom_insert.write_to_database(db_session);
	}

	return 0;
}

void
ResidueConformationFeatures::delete_record(
	StructureID struct_id,
	sessionOP db_session
){

	statement conformation_stmt(basic::database::safely_prepare_statement("DELETE FROM nonprotein_residue_conformation WHERE struct_id = ?;\n",db_session));
	conformation_stmt.bind(1,struct_id);
	basic::database::safely_write_to_database(conformation_stmt);
	statement angle_stmt(basic::database::safely_prepare_statement("DELETE FROM nonprotein_residue_angles WHERE struct_id = ?;\n",db_session));
	angle_stmt.bind(1,struct_id);
	basic::database::safely_write_to_database(angle_stmt);
	statement coords_stmt(basic::database::safely_prepare_statement("DELETE FROM residue_atom_coords WHERE struct_id = ?;",db_session));
	coords_stmt.bind(1,struct_id);
	basic::database::safely_write_to_database(coords_stmt);
	statement compact_coords_stmt(basic::database::safely_prepare_statement("DELETE FROM compact_residue_atom_coords WHERE struct_id = ?;",db_session));
	compact_coords_stmt.bind(1,struct_id);
	basic::database::safely_write_to_database(compact_coords_stmt);
}

void
ResidueConformationFeatures::load_into_pose(
	sessionOP db_session,
	StructureID struct_id,
	Pose & pose
){
	load_conformation(db_session, struct_id, pose);
}

void
ResidueConformationFeatures::load_conformation(
	sessionOP db_session,
	StructureID struct_id,
	Pose & pose
){
	if(!table_exists(db_session, "nonprotein_residue_conformation")) return;
	if(!table_exists(db_session, "nonprotein_residue_angles")) return;

	if(pose.is_fullatom()){

		std::string protein_string =
			"SELECT\n"
			"	seqpos,\n"
			"	phi,\n"
			"	psi,\n"
			"	omega\n"
			"FROM\n"
			"	nonprotein_residue_conformation\n"
			"WHERE\n"
			"	nonprotein_residue_conformation.struct_id=?;";

		statement protein_stmt(basic::database::safely_prepare_statement(protein_string,db_session));
		protein_stmt.bind(1,struct_id);

		std::string conformation_string =
			"SELECT\n"
			"	seqpos,\n"
			"	chinum,\n"
			"	chiangle\n"
			"FROM\n"
			"	nonprotein_residue_angles\n"
			"WHERE\n"
			"	nonprotein_residue_angles.struct_id=?;";

		statement conformation_stmt(basic::database::safely_prepare_statement(conformation_string,db_session));
		conformation_stmt.bind(1,struct_id);


		result protein_res(basic::database::safely_read_from_database(protein_stmt));
		while(protein_res.next()){
			Size seqpos;
			Real phi,psi,omega;
			protein_res >> seqpos >> phi >> psi >> omega;
			if (pose.residue_type(seqpos).is_protein()){
				pose.set_phi(seqpos,phi);
				pose.set_psi(seqpos,psi);
				pose.set_omega(seqpos,omega);
			}
		}
		result res_conformation(basic::database::safely_read_from_database(conformation_stmt));
		while(res_conformation.next()){
			//Size nchi(pose.residue_type(seqpos).nchi());
			Size seqpos;
			Size chinum;
			Real chiangle;
			if(compact_residue_schema_)
			{
				set_coords_for_residue_from_compact_schema(db_session,struct_id,seqpos,pose);
			}else
			{
				set_coords_for_residue(db_session,struct_id,seqpos,pose);

			}
			res_conformation >> seqpos >> chinum >> chiangle;
			pose.set_chi(chinum,seqpos,chiangle);
		}

	}else{

		std::string protein_string  =
			"SELECT\n"
			"	seqpos,\n"
			"	phi,\n"
			"	psi,\n"
			"	omega\n"
			"FROM\n"
			"	nonprotein_residue_conformation\n"
			"WHERE\n"
			"	nonprotein_residue_conformation.struct_id=?;";

		statement protein_stmt(basic::database::safely_prepare_statement(protein_string,db_session));
		protein_stmt.bind(1,struct_id);

		result protein_res(basic::database::safely_read_from_database(protein_stmt));
		while(protein_res.next()){
			Size seqpos;
			Real phi,psi,omega;
			protein_res >> seqpos >> phi >> psi >> omega;
			if (!pose.residue_type(seqpos).is_protein()){
				// WARNING why are you storing non-protein in the ProteinSilentReport?
				continue;
			}
			if(compact_residue_schema_)
			{
				set_coords_for_residue_from_compact_schema(db_session,struct_id,seqpos,pose);
			}else
			{
				set_coords_for_residue(db_session,struct_id,seqpos,pose);

			}
			pose.set_phi(seqpos,phi);
			pose.set_psi(seqpos,psi);
			pose.set_omega(seqpos,omega);
		}

	}
}

//This should be factored out into a non-member function.
void ResidueConformationFeatures::set_coords_for_residue(
		sessionOP db_session,
		StructureID struct_id,
		Size seqpos,
		Pose & pose
){

	std::string statement_string =
		"SELECT\n"
		"	atomno,\n"
		"	x,\n"
		"	y,\n"
		"	z\n"
		"FROM\n"
		"	residue_atom_coords\n"
		"WHERE\n"
		"residue_atom_coords.struct_id=? AND residue_atom_coords.seqpos=?;";

	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);
	stmt.bind(2,seqpos);

	result res(basic::database::safely_read_from_database(stmt));
	while(res.next()){
		Size atomno;
		Real x,y,z;
		res >> atomno >> x >> y >> z;

		core::id::AtomID atom_id(atomno,seqpos);
		core::Vector coords(x,y,z);
		pose.set_xyz(atom_id,coords);
	}


}

void
ResidueConformationFeatures::set_coords_for_residue_from_compact_schema(
	utility::sql_database::sessionOP db_session,
	StructureID struct_id,
	core::Size seqpos,
	core::pose::Pose & pose
) {

	std::string statement_string =
		"SELECT\n"
		"	coord_data,\n"
		"	atom_count\n"
		"FROM\n"
		"	compact_residue_atom_coords\n"
		"WHERE\n"
		"	compact_residue_atom_coords.struct_id=? AND compact_residue_atom_coords.seqpos=?";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);
	stmt.bind(2,seqpos);

	result res(basic::database::safely_read_from_database(stmt));
	while(res.next()){
		std::string coords;
		core::Size atom_count;
		res >> coords >>atom_count;


		utility::vector1<numeric::xyzVector<core::Real> > residue_coords(deserialize_xyz_coords(coords,atom_count)	);
		for(core::Size atomno = 1; atomno <= residue_coords.size();++atomno)
		{
			core::id::AtomID atom_id(atomno,seqpos);
			pose.set_xyz(atom_id,residue_coords[atomno]);
		}
	}
}

} // features
} // protocols
