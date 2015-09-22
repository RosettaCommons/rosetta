// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ProteinResidueConformationFeatures.cc
/// @brief  report idealized torsional DOFs Statistics Scientific Benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/ProteinResidueConformationFeatures.hh>

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
#include <utility/string_util.hh>
#include <utility/tools/make_vector.hh>
#include <protocols/features/util.hh>

//Basic Headers
#include <basic/database/sql_utils.hh>
#include <basic/Tracer.hh>
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

namespace protocols {
namespace features {

static THREAD_LOCAL basic::Tracer TR( "protocols.features.ProteinResidueConformationFeatures" );

using std::string;
using core::Size;
using core::Real;
using core::Vector;
using core::chemical::num_canonical_aas;
using core::conformation::Residue;
using core::pose::Pose;
using core::id::AtomID;
using utility::sql_database::sessionOP;
using utility::vector1;
using cppdb::statement;
using cppdb::result;
using basic::database::insert_statement_generator::InsertGenerator;
using basic::database::insert_statement_generator::RowDataBaseOP;
using basic::database::insert_statement_generator::RowData;


ProteinResidueConformationFeatures::ProteinResidueConformationFeatures()
{
	compact_residue_schema_ = basic::options::option[basic::options::OptionKeys::inout::dbms::use_compact_residue_schema]();
}

ProteinResidueConformationFeatures::ProteinResidueConformationFeatures(
	ProteinResidueConformationFeatures const & ) : FeaturesReporter()
{
	compact_residue_schema_ = basic::options::option[basic::options::OptionKeys::inout::dbms::use_compact_residue_schema]();
}


string
ProteinResidueConformationFeatures::type_name() const {
	return "ProteinResidueConformationFeatures";
}

void
ProteinResidueConformationFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const{

	using namespace basic::database::schema_generator;

	//******protein_residue_conformation******//
	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ), false);
	Column seqpos("seqpos", DbDataTypeOP( new DbInteger() ), false);
	Column secstruct("secstruct", DbDataTypeOP( new DbText() ), false);
	Column phi("phi", DbDataTypeOP( new DbDouble() ), false);
	Column psi("psi", DbDataTypeOP( new DbDouble() ), false);
	Column omega("omega", DbDataTypeOP( new DbDouble() ), false);
	Column chi1("chi1", DbDataTypeOP( new DbDouble() ), false);
	Column chi2("chi2", DbDataTypeOP( new DbDouble() ), false);
	Column chi3("chi3", DbDataTypeOP( new DbDouble() ), false);
	Column chi4("chi4", DbDataTypeOP( new DbDouble() ), false);


	utility::vector1<Column> prot_res_pkeys;
	prot_res_pkeys.push_back(struct_id);
	prot_res_pkeys.push_back(seqpos);

	utility::vector1<Column> fkey_cols;
	fkey_cols.push_back(struct_id);
	fkey_cols.push_back(seqpos);

	utility::vector1<std::string> fkey_reference_cols;
	fkey_reference_cols.push_back("struct_id");
	fkey_reference_cols.push_back("resNum");

	Schema protein_residue_conformation("protein_residue_conformation", PrimaryKey(prot_res_pkeys));
	protein_residue_conformation.add_column(struct_id);
	protein_residue_conformation.add_column(seqpos);
	protein_residue_conformation.add_column(secstruct);
	protein_residue_conformation.add_column(phi);
	protein_residue_conformation.add_column(psi);
	protein_residue_conformation.add_column(omega);
	protein_residue_conformation.add_column(chi1);
	protein_residue_conformation.add_column(chi2);
	protein_residue_conformation.add_column(chi3);
	protein_residue_conformation.add_column(chi4);
	protein_residue_conformation.add_foreign_key(ForeignKey(fkey_cols, "residues", fkey_reference_cols, true));

	protein_residue_conformation.write(db_session);

	if ( compact_residue_schema_ ) {
		//******compact_residue_atom_coords*****//
		Column coord_data("coord_data",DbDataTypeOP( new DbText() ));
		Column atom_count("atom_count",DbDataTypeOP( new DbInteger() ),false);
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

	} else {
		//******residue_atom_coords******//
		Column atomno("atomno", DbDataTypeOP( new DbInteger() ), false);
		Column x("x", DbDataTypeOP( new DbDouble() ), false);
		Column y("y", DbDataTypeOP( new DbDouble() ), false);
		Column z("z", DbDataTypeOP( new DbDouble() ), false);

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
ProteinResidueConformationFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	return dependencies;
}


Size
ProteinResidueConformationFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID struct_id,
	sessionOP db_session
){
	bool fullatom(pose.is_fullatom());

	// Check to see if this structure is ideal
	// KAB - this code has also been copied to PoseConformationFeatures
	// so that that features reporter can save if the pose is ideal.
	// Theoretically, it could be possible to remove the duplication,
	// but it would probably be have to cached in pose or somewhere,
	// and as force_nonideal_structure is the default option, it is unlikely
	// this code will not be run very often
	bool ideal = true;
	if ( !basic::options::option[basic::options::OptionKeys::out::file::force_nonideal_structure]() ) {
		core::conformation::Conformation const & conformation(pose.conformation());
		for ( core::Size resn=1; resn <= pose.n_residue(); ++resn ) {
			if ( !check_relevant_residues( relevant_residues, resn ) ) continue;
			bool residue_status(core::conformation::is_ideal_position(resn,conformation));
			if ( !residue_status ) {
				ideal = false;
				break;
			}
		}
	} else {
		ideal = false;
	}

	InsertGenerator conformation_insert("protein_residue_conformation");
	conformation_insert.add_column("struct_id");
	conformation_insert.add_column("seqpos");
	conformation_insert.add_column("secstruct");
	conformation_insert.add_column("phi");
	conformation_insert.add_column("psi");
	conformation_insert.add_column("omega");
	conformation_insert.add_column("chi1");
	conformation_insert.add_column("chi2");
	conformation_insert.add_column("chi3");
	conformation_insert.add_column("chi4");

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

	RowDataBaseOP struct_id_data( new RowData<StructureID>("struct_id",struct_id) );

	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( !check_relevant_residues( relevant_residues, i ) ) continue;
		Residue const & resi = pose.residue(i);
		if ( resi.aa() > num_canonical_aas ) {
			continue;
		}
		std::string secstruct = utility::to_string<char>(pose.secstruct(i));

		Real phi = resi.mainchain_torsion(1);
		Real psi = resi.mainchain_torsion(2);
		Real omega = resi.mainchain_torsion(3);

		Real chi1 = fullatom && resi.nchi() >= 1 ? resi.chi(1) : 0.0;
		Real chi2 = fullatom && resi.nchi() >= 2 ? resi.chi(2) : 0.0;
		Real chi3 = fullatom && resi.nchi() >= 3 ? resi.chi(3) : 0.0;
		Real chi4 = fullatom && resi.nchi() >= 4 ? resi.chi(4) : 0.0;

		RowDataBaseOP seqpos_data( new RowData<Size>("seqpos",i) );
		RowDataBaseOP secstruct_data( new RowData<std::string>("secstruct",secstruct) );
		RowDataBaseOP phi_data( new RowData<Real>("phi",phi) );
		RowDataBaseOP psi_data( new RowData<Real>("psi",psi) );
		RowDataBaseOP omega_data( new RowData<Real>("omega",omega) );
		RowDataBaseOP chi1_data( new RowData<Real>("chi1",chi1) );
		RowDataBaseOP chi2_data( new RowData<Real>("chi2",chi2) );
		RowDataBaseOP chi3_data( new RowData<Real>("chi3",chi3) );
		RowDataBaseOP chi4_data( new RowData<Real>("chi4",chi4) );


		conformation_insert.add_row(utility::tools::make_vector(
			struct_id_data,seqpos_data,secstruct_data,phi_data,
			psi_data,omega_data,chi1_data,chi2_data,chi3_data,chi4_data));

		// KAB - We're still only saving coordinate data if the structure is not ideal
		// If anyone is interested in saving this data always, you probably can make it happen now that we save if a structure
		// is ideal or not in the pose_conformations table. This is what I did in order to always save backbone torsions.
		if ( !ideal ) {
			if ( compact_residue_schema_ ) {
				std::string residue_data_string(serialize_residue_xyz_coords(resi));
				RowDataBaseOP atom_count( new RowData<core::Size>("atom_count",resi.natoms()) );
				RowDataBaseOP residue_data( new RowData<std::string>("coord_data",residue_data_string) );
				compact_residue_insert.add_row(utility::tools::make_vector(struct_id_data,atom_count,seqpos_data,residue_data));
			} else {
				for ( Size atom = 1; atom <= resi.natoms(); ++atom ) {
					core::Vector coords = resi.xyz(atom);

					RowDataBaseOP atom_data( new RowData<Size>("atomno",atom) );
					RowDataBaseOP x_data( new RowData<Real>("x",coords.x()) );
					RowDataBaseOP y_data( new RowData<Real>("y",coords.y()) );
					RowDataBaseOP z_data( new RowData<Real>("z",coords.z()) );

					atom_insert.add_row(utility::tools::make_vector(
						struct_id_data,seqpos_data,atom_data,x_data,y_data,z_data));
				}
			}
		}
	}

	conformation_insert.write_to_database(db_session);
	if ( compact_residue_schema_ ) {
		compact_residue_insert.write_to_database(db_session);
	} else {
		atom_insert.write_to_database(db_session);
	}

	return 0;
}

void
ProteinResidueConformationFeatures::delete_record(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session
){

	//std::string struct_id_string(to_string(struct_id));
	statement conf_stmt(basic::database::safely_prepare_statement("DELETE FROM protein_residue_conformation WHERE struct_id = ?;\n",db_session));
	conf_stmt.bind(1,struct_id);
	basic::database::safely_write_to_database(conf_stmt);
	if ( compact_residue_schema_ ) {
		statement compact_coords_stmt(basic::database::safely_prepare_statement("DELETE FROM compact_residue_atom_coords WHERE struct_id = ?;",db_session));
		compact_coords_stmt.bind(1,struct_id);
		basic::database::safely_write_to_database(compact_coords_stmt);
	} else {
		statement atom_stmt(basic::database::safely_prepare_statement("DELETE FROM residue_atom_coords WHERE struct_id = ?;\n",db_session));
		atom_stmt.bind(1,struct_id);
		basic::database::safely_write_to_database(atom_stmt);
	}
}


void
ProteinResidueConformationFeatures::load_into_pose(
	sessionOP db_session,
	StructureID struct_id,
	Pose & pose,
	bool ideal
){
	load_conformation(db_session, struct_id, pose, ideal);
}

void
ProteinResidueConformationFeatures::load_conformation(
	sessionOP db_session,
	StructureID struct_id,
	Pose & pose,
	bool ideal
){

	if ( !basic::database::table_exists(db_session, "protein_residue_conformation") ) {
		TR << "WARNING: protein_residue_conformation table does not exist and thus respective data will not be added to the pose!" << std::endl;
		return;
	}

	if ( compact_residue_schema_ ) {
		set_coords_for_residue_from_compact_schema(db_session,struct_id,pose);
	} else {
		set_coords_for_residues(db_session,struct_id,pose);
	}


	if ( pose.is_fullatom() ) {
		std::string statement_string =
			"SELECT\n"
			"\tsecstruct,\n"
			"\tphi,\n"
			"\tpsi,\n"
			"\tomega,\n"
			"\tchi1,\n"
			"\tchi2,\n"
			"\tchi3,\n"
			"\tchi4\n"
			"FROM\n"
			"\tprotein_residue_conformation\n"
			"WHERE\n"
			"\tprotein_residue_conformation.struct_id=?;";
		statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
		//std::string struct_id_string(to_string(struct_id));
		stmt.bind(1,struct_id);

		result res(basic::database::safely_read_from_database(stmt));

		Size seqpos(1);
		while ( res.next() ) {
			Real phi,psi,omega,chi1,chi2,chi3,chi4;
			std::string secstruct;
			res >> secstruct >> phi >> psi >> omega >> chi1 >> chi2 >> chi3 >> chi4 ;
			if ( !pose.residue_type(seqpos).is_protein() ) {
				// WARNING why are you storing non-protein in the ProteinSilentReport?
				continue;
			}
			if ( ideal && !(phi < 0.00001 && psi < 0.00001 && omega < 0.00001) ) {
				pose.set_phi(seqpos,phi);
				pose.set_psi(seqpos,psi);
				pose.set_omega(seqpos,omega);
			}
			pose.set_secstruct(seqpos,static_cast<char>(secstruct[0]));
			Size nchi(pose.residue_type(seqpos).nchi());
			if ( 1 <= nchi ) pose.set_chi(1, seqpos, chi1);
			if ( 2 <= nchi ) pose.set_chi(2, seqpos, chi2);
			if ( 3 <= nchi ) pose.set_chi(3, seqpos, chi3);
			if ( 4 <= nchi ) pose.set_chi(4, seqpos, chi4);
			seqpos++;
		}

	} else {
		// statement stmt = (*db_session) <<
		//  "SELECT\n"
		//  " seqpos,\n"
		//  " phi,\n"
		//  " psi,\n"
		//  " omega\n"
		//  "FROM\n"
		//  " protein_residue_conformation\n"
		//  "WHERE\n"
		//  " protein_residue_conformation.struct_id=?;" << struct_id;
		//
		// result res(basic::database::safely_read_from_database(stmt));
		// while(res.next()){
		//  Size seqpos;
		//  Real phi,psi,omega;
		//  res >> seqpos >> phi >> psi >> omega;
		//  if (!pose.residue_type(seqpos).is_protein()){
		//   // WARNING why are you storing non-protein in the ProteinSilentReport?
		//   continue;
		//  }
		//
		//  //pose.set_phi(seqpos,phi);
		//  //pose.set_psi(seqpos,psi);
		//  //pose.set_omega(seqpos,omega);
		// }
	}
}

void
ProteinResidueConformationFeatures::check_num_requested_atoms(
	Size num_requested_atoms,
	Size pose_resNum,
	Pose & pose,
	Size resNum,
	StructureID struct_id,
	sessionOP db_session
) const {
	if ( num_requested_atoms != pose.residue_type(pose_resNum).natoms() ) {
		string get_res_type_stmt_str("SELECT res_type FROM residues WHERE struct_id = ? AND resNum = ?;");
		statement get_res_type_stmt(basic::database::safely_prepare_statement(get_res_type_stmt_str, db_session));
		get_res_type_stmt.bind(1, struct_id);
		get_res_type_stmt.bind(2, resNum);
		result get_res_type_res(basic::database::safely_read_from_database(get_res_type_stmt));
		get_res_type_res.next();
		string db_res_type;
		get_res_type_res >> db_res_type;

		if (
				db_res_type == "CYD" && num_requested_atoms == 10 &&
				pose.residue_type(pose_resNum).name() == "CYS"
				) {
			// The pose incorrectly says this cystine is
			// forming a disulfide while the coordinates say that is
			// not. Convert the residue type to CYD.
			bool success = core::conformation::change_cys_state(
				pose_resNum, "CYD", pose.conformation());
			if ( !success ) {
				std::stringstream errmsg;
				errmsg
					<< "Failed to convert cystine from CYS to CYD for residue '" << resNum << "'";
				if ( pose_resNum != resNum ) {
					errmsg
						<< ", which has been renumbered to '" << pose_resNum << "'." << std::endl;
				} else {
					errmsg
						<< "." << std::endl;
				}
				utility_exit_with_message(errmsg.str());
			} else {
				// should probably rebuild disulfides when done...
				return;
			}
		}


		std::stringstream errmsg;
		errmsg
			<< "Attempting to set atom coordinates for residue '" << resNum << "'";
		if ( pose_resNum != resNum ) {
			errmsg
				<< ", which has been renumbered to '" << pose_resNum << "'." << std::endl;
		} else {
			errmsg
				<< "." << std::endl;
		}

		errmsg
			<< "The pose residue has '" << pose.residue_type(pose_resNum).natoms() << "' coordinates but '" << num_requested_atoms << "' were provided." << std::endl;
		if ( db_res_type != pose.residue_type(pose_resNum).name() ) {
			errmsg
				<< "The pose residue type is '" << pose.residue_type(pose_resNum).name() << "' while the residue type in the database is '" << db_res_type << "'." << std::endl;
		} else {
			errmsg
				<< "The residue type in the database is '" << db_res_type << "' and it matches the pose residue type." << std::endl;
		}
		utility_exit_with_message(errmsg.str());
	}
}

//This should be factored out into a non-member function.
/// This is a little complicated:
/// PoseConformationFeatures: renumbers the pose then reports the relevant residues
/// ResidueFeatures: doesn't renumber the pose and reports the relevant residues
/// ProteinResidueConformationFeatures: don't renumber the pose but also only report features for (canonical) protein residues
///
/// This means that to reconstruct the pose conformation, first the conformations
/// stored in residue_atom_coords must be aligned with residues table
/// to get the appropriate gaps for non-canonical residues, then these
/// need to be renumbered to align with the renumbered pose reported
/// in PoseConformationFEatures
///
///    pose_resNum -> the numbering in the pose to be filled
///    resNum -> the numbering in the residues table

void ProteinResidueConformationFeatures::set_coords_for_residues(
	sessionOP db_session,
	StructureID struct_id,
	Pose & pose
){
	// lookup and set all the atoms at once because each query is
	// roughly O(n*log(n) + k) where n is the size of the tables and k
	// is the number of rows returned. Doing it all at once means you
	// only have to pay the n*log(n) cost once.

	std::string statement_string =
		"SELECT\n"
		"\tr.resNum,\n"
		"\tc.atomno,\n"
		"\tc.x,\n"
		"\tc.y,\n"
		"\tc.z\n"
		"FROM\n"
		"\tresidues AS r LEFT JOIN\n"
		"\tresidue_atom_coords AS c ON\n"
		"\tr.struct_id = c.struct_id AND r.resNum = c.seqpos\n"
		"WHERE\n"
		"\tr.struct_id=?;";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);

	result res(basic::database::safely_read_from_database(stmt));

	vector1< AtomID > atom_ids;
	vector1< Vector > coords;
	vector1< Size > missing_coordinates;
	bool initial_iteration(true);
	bool last_residue_missing(true);
	Size pose_resNum(1);
	Size prev_resNum(0);
	Size natoms_in_residue(0);
	Size resNum, atomno;
	Real x, y, z;
	while ( res.next() ) {
		res >> resNum >> atomno >> x >> y >> z;
		if ( !initial_iteration && prev_resNum != resNum ) {
			if ( !last_residue_missing ) {
				check_num_requested_atoms(natoms_in_residue, pose_resNum, pose, prev_resNum, struct_id, db_session);
			}
			pose_resNum++;
			natoms_in_residue = 0;
		}

		if ( res.is_null(1) ) {
			missing_coordinates.push_back(pose_resNum);
			last_residue_missing = true;
			continue;
		} else {
			last_residue_missing = false;
		}
		atom_ids.push_back(AtomID(atomno, pose_resNum));
		coords.push_back(Vector(x,y,z));
		natoms_in_residue++;
		prev_resNum = resNum;
		initial_iteration=false;
	}
	if ( !last_residue_missing ) {
		check_num_requested_atoms(natoms_in_residue, pose_resNum, pose, resNum, struct_id, db_session);
	}

	// use the batch_set_xyz because it doesn't trigger a coordinate
	// update after setting each atom.
	pose.batch_set_xyz(atom_ids, coords);


	if ( missing_coordinates.size() > 0 ) {
		TR.Warning << "In loading the residue coodinates, some of the residues did not have coordinates specified:" << std::endl << "\t[";
		for ( Size ii=1; ii <= missing_coordinates.size(); ++ii ) {
			TR.Warning << missing_coordinates[ii] << ",";
		}
		TR.Warning << "]" << std::endl;
		TR.Warning << "This can happen because you are using ProteinResidueConformationFeatures with a structure that contains non-protein residues, for example. To avoid this, extract with the ResidueConformationFeatures instead." << std::endl;
		//  TR.Warning << "These residues will be deleted from the pose." << std::endl;
		//  for(Size ii=1; ii < missing_coordinates.size(); ++ii){
		//   pose.conformation().delete_residue_slow(ii);
		//  }
	}
}

void
ProteinResidueConformationFeatures::set_coords_for_residue_from_compact_schema(
	utility::sql_database::sessionOP db_session,
	StructureID struct_id,
	core::pose::Pose & pose
) {

	std::string statement_string =
		"SELECT\n"
		"\tcoord_data,\n"
		"\tatom_count\n"
		"FROM\n"
		"\tcompact_residue_atom_coords\n"
		"WHERE\n"
		"\tcompact_residue_atom_coords.struct_id=?";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);

	result res(basic::database::safely_read_from_database(stmt));
	Size seqpos(1);
	while ( res.next() ) {
		std::string coords;
		core::Size atom_count = 0;
		res >> coords >> atom_count;
		utility::vector1<numeric::xyzVector<core::Real> > residue_coords(deserialize_xyz_coords(coords,atom_count) );
		for ( core::Size atomno = 1; atomno <= residue_coords.size(); ++atomno ) {
			core::id::AtomID atom_id(atomno,seqpos);
			pose.set_xyz(atom_id,residue_coords[atomno]);
		}
		seqpos++;
	}
}

} // namespace
} // namespace
