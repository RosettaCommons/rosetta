// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/PoseConformationFeatures.cc
/// @brief  report comments stored with each pose
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/PoseConformationFeatures.hh>

//External

// Project Headers
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/io/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/types.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/Constraint.hh>
#include <basic/database/schema_generator/DbDataType.hh>

#include <basic/database/insert_statement_generator/InsertGenerator.hh>
#include <basic/database/insert_statement_generator/RowData.hh>

// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/database/sql_utils.hh>
#include <basic/Tracer.hh>


// Numeric Headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/tools/make_vector.hh>


// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <cmath>
#include <sstream>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/features/feature_schemas.hh>
#include <protocols/features/PoseConformationFeaturesCreator.hh>

namespace protocols {
namespace features {

static THREAD_LOCAL basic::Tracer TR( "protocols.features.PoseConformationFeatures" );

using std::string;
using std::stringstream;
using std::endl;
using core::Real;
using core::Size;
using core::Vector;
using core::kinematics::FoldTree;
using core::kinematics::Jump;
using core::pose::Pose;
using core::pose::PoseOP;
using core::pose::PoseCOP;
using core::pose::make_pose_from_sequence;
using core::chemical::ResidueTypeSetCOP;
using core::chemical::ResidueTypeCOPs;
using core::chemical::ChemicalManager;
using core::chemical::FA_STANDARD;
using core::chemical::CENTROID;
using core::io::pose_from_pose;
using core::kinematics::FoldTree;
using core::kinematics::Edge;
using core::kinematics::Jump;
using core::kinematics::RT;
using numeric::xyzMatrix;
using numeric::xyzVector;
using utility::vector1;
using utility::sql_database::DatabaseSessionManager;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;
using basic::database::insert_statement_generator::InsertGenerator;
using basic::database::insert_statement_generator::RowDataBaseOP;
using basic::database::insert_statement_generator::RowData;

// XRW TEMP string
// XRW TEMP PoseConformationFeatures::type_name() const { return "PoseConformationFeatures"; }

void
PoseConformationFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const{

	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ), false /*not null*/, false /*don't autoincrement*/);

	/******pose_conformations******/
	Schema pose_conformations("pose_conformations");
	pose_conformations.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));

	pose_conformations.add_column( Column("annotated_sequence", DbDataTypeOP( new DbText() )) );
	pose_conformations.add_column( Column("total_residue", DbDataTypeOP( new DbInteger() )) );
	pose_conformations.add_column( Column("fullatom", DbDataTypeOP( new DbInteger() )) );
	pose_conformations.add_column( Column("ideal", DbDataTypeOP( new DbInteger() )) ); // 0 if not ideal, 1 if ideal

	pose_conformations.write(db_session);

	/******fold_trees******/
	Schema fold_trees("fold_trees");
	fold_trees.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));

	fold_trees.add_column( Column("start_res", DbDataTypeOP( new DbInteger() )) );
	fold_trees.add_column( Column("start_atom", DbDataTypeOP( new DbText() )) );
	fold_trees.add_column( Column("stop_res", DbDataTypeOP( new DbInteger() )) );
	fold_trees.add_column( Column("stop_atom", DbDataTypeOP( new DbText() )) );
	fold_trees.add_column( Column("label", DbDataTypeOP( new DbInteger() )) );
	fold_trees.add_column( Column("keep_stub_in_residue", DbDataTypeOP( new DbInteger() )) );

	fold_trees.write(db_session);

	/******jumps******/
	Schema jumps("jumps");
	jumps.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));
	jumps.add_column( Column("jump_id", DbDataTypeOP( new DbInteger() )) );
	jumps.add_column( Column("xx", DbDataTypeOP( new DbDouble() )) );
	jumps.add_column( Column("xy", DbDataTypeOP( new DbDouble() )) );
	jumps.add_column( Column("xz", DbDataTypeOP( new DbDouble() )) );
	jumps.add_column( Column("yx", DbDataTypeOP( new DbDouble() )) );
	jumps.add_column( Column("yy", DbDataTypeOP( new DbDouble() )) );
	jumps.add_column( Column("yz", DbDataTypeOP( new DbDouble() )) );
	jumps.add_column( Column("zx", DbDataTypeOP( new DbDouble() )) );
	jumps.add_column( Column("zy", DbDataTypeOP( new DbDouble() )) );
	jumps.add_column( Column("zz", DbDataTypeOP( new DbDouble() )) );
	jumps.add_column( Column("x", DbDataTypeOP( new DbDouble() )) );
	jumps.add_column( Column("y", DbDataTypeOP( new DbDouble() )) );
	jumps.add_column( Column("z", DbDataTypeOP( new DbDouble() )) );

	jumps.write(db_session);

	/******chain_endings******/
	Schema chain_endings("chain_endings");
	chain_endings.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));
	chain_endings.add_column( Column("end_pos", DbDataTypeOP( new DbInteger() )) );

	chain_endings.write(db_session);
}

utility::vector1<std::string>
PoseConformationFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("StructureFeatures");
	return dependencies;
}

Size
PoseConformationFeatures::report_features(
	Pose const & pose_orig,
	vector1< bool > const & relevant_residues,
	StructureID struct_id,
	sessionOP db_session
){
	vector1< Size > residue_indices;
	for ( Size i = 1; i <= relevant_residues.size(); ++i ) {
		if ( relevant_residues[i] ) residue_indices.push_back(i);
	}

	if ( residue_indices.size() == pose_orig.size() ) {
		return report_features_implementation(pose_orig, struct_id, db_session);
	}
	// else...
	Pose pose;
	pose_from_pose( pose, pose_orig, residue_indices );
	return report_features_implementation(pose, struct_id, db_session);
}

Size
PoseConformationFeatures::report_features_implementation(
	Pose const & pose,
	StructureID struct_id,
	sessionOP db_session
){
	FoldTree const & fold_tree(pose.conformation().fold_tree());
	//assume non-trivial fold_tree only if more than one edge, i.e., EDGE 1 <nres> -1

	InsertGenerator fold_tree_insert("fold_trees");
	fold_tree_insert.add_column("struct_id");
	fold_tree_insert.add_column("start_res");
	fold_tree_insert.add_column("start_atom");
	fold_tree_insert.add_column("stop_res");
	fold_tree_insert.add_column("stop_atom");
	fold_tree_insert.add_column("label");
	fold_tree_insert.add_column("keep_stub_in_residue");

	RowDataBaseOP struct_id_data( new RowData<StructureID>("struct_id",struct_id) );

	for ( auto const & it : fold_tree ) {

		int start_res(it.start()), stop_res(it.stop()), label(it.label());
		string start_atom(it.start_atom()), stop_atom(it.stop_atom());
		bool keep_stub_in_residue(it.keep_stub_in_residue());


		RowDataBaseOP start_res_data( new RowData<int>("start_res",start_res) );
		RowDataBaseOP start_atom_data( new RowData<string>("start_atom",start_atom) );
		RowDataBaseOP stop_res_data( new RowData<int>("stop_res",stop_res) );
		RowDataBaseOP stop_atom_data( new RowData<string>("stop_atom",stop_atom) );
		RowDataBaseOP label_data( new RowData<int>("label",label) );
		RowDataBaseOP keep_stub_data( new RowData<bool>("keep_stub_in_residue",keep_stub_in_residue) );

		fold_tree_insert.add_row(
			utility::tools::make_vector(struct_id_data,start_res_data,start_atom_data,stop_res_data,stop_atom_data,label_data,keep_stub_data));
	}

	fold_tree_insert.write_to_database(db_session);

	InsertGenerator jump_insert("jumps");
	jump_insert.add_column("struct_id");
	jump_insert.add_column("jump_id");
	jump_insert.add_column("xx");
	jump_insert.add_column("xy");
	jump_insert.add_column("xz");
	jump_insert.add_column("yx");
	jump_insert.add_column("yy");
	jump_insert.add_column("yz");
	jump_insert.add_column("zx");
	jump_insert.add_column("zy");
	jump_insert.add_column("zz");
	jump_insert.add_column("x");
	jump_insert.add_column("y");
	jump_insert.add_column("z");

	for ( Size nr = 1; nr <= fold_tree.num_jump(); nr++ ) {
		Jump const & jump(pose.jump(nr));
		xyzMatrix< Real > const & r(jump.get_rotation());
		Real xx(r.xx()), xy(r.xy()), xz(r.xz());
		Real yx(r.yx()), yy(r.yy()), yz(r.yz());
		Real zx(r.zx()), zy(r.zy()), zz(r.zz());
		Vector const & t(jump.get_translation());
		Real x(t.x()), y(t.y()), z(t.z());

		RowDataBaseOP jump_id_data( new RowData<Size>("jump_id",nr) );
		RowDataBaseOP xx_data( new RowData<Real>("xx",xx) );
		RowDataBaseOP xy_data( new RowData<Real>("xy",xy) );
		RowDataBaseOP xz_data( new RowData<Real>("xz",xz) );
		RowDataBaseOP yx_data( new RowData<Real>("yx",yx) );
		RowDataBaseOP yy_data( new RowData<Real>("yy",yy) );
		RowDataBaseOP yz_data( new RowData<Real>("yz",yz) );
		RowDataBaseOP zx_data( new RowData<Real>("zx",zx) );
		RowDataBaseOP zy_data( new RowData<Real>("zy",zy) );
		RowDataBaseOP zz_data( new RowData<Real>("zz",zz) );
		RowDataBaseOP x_data( new RowData<Real>("x",x) );
		RowDataBaseOP y_data( new RowData<Real>("y",y) );
		RowDataBaseOP z_data( new RowData<Real>("z",z) );
		jump_insert.add_row(
			utility::tools::make_vector(struct_id_data,jump_id_data,xx_data,xy_data,xz_data,yx_data,yy_data,yz_data,zx_data,zy_data,zz_data,x_data,y_data,z_data));
	}
	jump_insert.write_to_database(db_session);

	InsertGenerator chain_ending_insert("chain_endings");
	chain_ending_insert.add_column("struct_id");
	chain_ending_insert.add_column("end_pos");
	for ( Size const end_pos : pose.conformation().chain_endings() ) {
		RowDataBaseOP end_pos_data( new RowData<Size>("end_pos",end_pos) );

		chain_ending_insert.add_row(
			utility::tools::make_vector(struct_id_data,end_pos_data));
	}
	chain_ending_insert.write_to_database(db_session);

	string annotated_sequence(pose.annotated_sequence(true));


	InsertGenerator pose_conformation_insert("pose_conformations");
	pose_conformation_insert.add_column("struct_id");
	pose_conformation_insert.add_column("annotated_sequence");
	pose_conformation_insert.add_column("total_residue");
	pose_conformation_insert.add_column("fullatom");
	pose_conformation_insert.add_column("ideal");

	RowDataBaseOP annotated_sequence_data( new RowData<string>("annotated_sequence",annotated_sequence) );
	RowDataBaseOP total_residue_data( new RowData<Size>("total_residue",pose.size()) );
	RowDataBaseOP fullatom_data( new RowData<bool>("fullatom",pose.is_fullatom()) );

	// KAB -
	// Determine if pose is ideal
	// Saving this here lets the database pose loader know later if it should load backbone torsions
	// This code and the comment below were copied from ProteinResidueConformationFeatures

	//If you have a non ideal structure, and you store both cartesian coordinates and backbone
	//chi angles and read them into a pose, most of the backbone oxygens will be placed in correctly
	//I currently have absolutely no idea why this is, but this fixes it.
	//It is worth noting that the current implementation of Binary protein silent files does the same thing
	bool ideal = true;
	if ( !basic::options::option[basic::options::OptionKeys::out::file::force_nonideal_structure]() ) {
		core::conformation::Conformation const & conformation(pose.conformation());
		for ( core::Size resn=1; resn <= pose.size(); ++resn ) {
			bool residue_status = core::conformation::is_ideal_position(resn,conformation);
			if ( !residue_status ) {
				ideal = false;
				break;
			}
		}
	} else {
		ideal = false;
	}
	RowDataBaseOP ideal_data( new RowData<bool>("ideal", ideal) );

	pose_conformation_insert.add_row(
		utility::tools::make_vector(struct_id_data,annotated_sequence_data,total_residue_data,fullatom_data, ideal_data));

	pose_conformation_insert.write_to_database(db_session);
	return 0;
}

void PoseConformationFeatures::delete_record(
	StructureID struct_id,
	sessionOP db_session
){

	statement conf_stmt(basic::database::safely_prepare_statement("DELETE FROM pose_conformations WHERE struct_id = ?;\n",db_session));
	conf_stmt.bind(1,struct_id);
	basic::database::safely_write_to_database(conf_stmt);

	statement fold_stmt(basic::database::safely_prepare_statement("DELETE FROM fold_trees WHERE struct_id = ?;\n",db_session));
	fold_stmt.bind(1,struct_id);
	basic::database::safely_write_to_database(fold_stmt);

	statement jump_stmt(basic::database::safely_prepare_statement("DELETE FROM jumps WHERE struct_id = ?;\n",db_session));
	jump_stmt.bind(1,struct_id);
	basic::database::safely_write_to_database(jump_stmt);

	statement chain_stmt(basic::database::safely_prepare_statement("DELETE FROM chain_endings WHERE struct_id = ?;",db_session));
	chain_stmt.bind(1,struct_id);
	basic::database::safely_write_to_database(chain_stmt);
}

void
PoseConformationFeatures::load_into_pose(
	sessionOP db_session,
	StructureID struct_id,
	Pose & pose,
	bool & ideal
){
	ideal = load_sequence(db_session, struct_id, pose);
	load_fold_tree(db_session, struct_id, pose);
	load_jumps(db_session, struct_id, pose);
	load_chain_endings(db_session, struct_id, pose);
}

bool
PoseConformationFeatures::load_sequence(
	sessionOP db_session,
	StructureID struct_id,
	Pose & pose
){

	if ( !basic::database::table_exists(db_session, "pose_conformations") ) {
		TR.Warning << "pose_conformations table does not exist and thus respective data will not be added to the pose!" << std::endl;
		return true; // Assume structure is ideal
	}


	string annotated_sequence;
	Size total_residue, fullatom, ideal;
	{
		std::string statement_string =
			"SELECT\n"
			"\tannotated_sequence,\n"
			"\ttotal_residue,\n"
			"\tfullatom,\n"
			" ideal\n"
			"FROM\n"
			"\tpose_conformations\n"
			"WHERE\n"
			"\tpose_conformations.struct_id = ?;";
		statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
		stmt.bind(1,struct_id);
		result res(basic::database::safely_read_from_database(stmt));

		if ( !res.next() ) {
			stringstream error_message;
			error_message << "Unable to locate structure with struct_id '" << struct_id << "'";
			utility_exit_with_message(error_message.str());
		}
		res >> annotated_sequence >> total_residue >> fullatom >> ideal;
	}
	ResidueTypeSetCOP residue_set(
		ChemicalManager::get_instance()->residue_type_set(
		fullatom ? FA_STANDARD : CENTROID));


	if ( basic::database::table_exists(db_session, "residues") ) {
		std::string statement_string =
			"SELECT\n"
			"\tres_type\n"
			"FROM\n"
			"\tresidues\n"
			"WHERE\n"
			"\tstruct_id = ?;";
		statement stmt(
			basic::database::safely_prepare_statement(statement_string,db_session));
		stmt.bind(1,struct_id);
		result res(basic::database::safely_read_from_database(stmt));

		ResidueTypeCOPs requested_types;
		while ( res.next() ) {
			string res_type;
			res >> res_type;

			if ( res_type=="CYD" ) {
				res_type="CYS:disulfide"; //JAB - Backwards compatability with older databases
			}
			requested_types.push_back( residue_set->name_mapOP( res_type ) );
		}
		make_pose_from_sequence(pose, requested_types);

	} else {

		make_pose_from_sequence(pose, annotated_sequence, *residue_set);
		runtime_assert(pose.size() == total_residue );
	}

	return ideal;
}


void
PoseConformationFeatures::load_fold_tree(
	sessionOP db_session,
	StructureID struct_id,
	Pose & pose
){

	if ( !basic::database::table_exists(db_session, "fold_trees") ) {
		TR.Warning << "fold_trees table does not exist and thus respective data will not be added to the pose!" << std::endl;
		return;
	}

	statement stmt = (*db_session) <<
		"SELECT\n"
		"\tstart_res,\n"
		"\tstart_atom,\n"
		"\tstop_res,\n"
		"\tstop_atom,\n"
		"\tlabel,\n"
		"\tkeep_stub_in_residue\n"
		"FROM\n"
		"\tfold_trees\n"
		"WHERE\n"
		"\tfold_trees.struct_id=?;" << struct_id;

	result res(basic::database::safely_read_from_database(stmt));

	FoldTree t = FoldTree();
	while ( res.next() ) {
		int start_res, stop_res, label;
		string start_atom, stop_atom;
		int keep_stub_in_residue;
		res >> start_res >> start_atom >> stop_res >> stop_atom >> label >> keep_stub_in_residue;
		if ( label == -2 || label > 0 ) { //CHEMICAL or JUMP
			t.add_edge(Edge(
				start_res, stop_res, label, start_atom, stop_atom, keep_stub_in_residue));
		} else {
			t.add_edge(Edge(start_res, stop_res, label, "", "", keep_stub_in_residue));
		}
	}
	// TODO verify that pose.fold_tree(t) is ok (not cleared from the stack)
	pose.fold_tree(t);

	TR.Debug << "Fold tree loaded" << std::endl;
}

void
PoseConformationFeatures::load_jumps(
	sessionOP db_session,
	StructureID struct_id,
	Pose & pose
){
	if ( !basic::database::table_exists(db_session, "jumps") ) {
		TR.Warning << "jumps table does not exist and thus respective data will not be added to the pose!" << std::endl;
		return;
	}

	//note the Conformation object sorts the chain_endings after they are passed in.
	std::string statement_string =
		"SELECT\n"
		"\tjump_id,\n"
		"\txx,\n"
		"\txy,\n"
		"\txz,\n"
		"\tyx,\n"
		"\tyy,\n"
		"\tyz,\n"
		"\tzx,\n"
		"\tzy,\n"
		"\tzz,\n"
		"\tx,\n"
		"\ty,\n"
		"\tz\n"
		"FROM\n"
		"\tjumps\n"
		"WHERE\n"
		"\tjumps.struct_id=?;";
	statement stmt(
		basic::database::safely_prepare_statement(statement_string,db_session));

	stmt.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(stmt));
	while ( res.next() ) {
		Size jump_id;
		Real xx, xy, xz, yx, yy, yz, zx, zy, zz, x, y, z;
		res >> jump_id;
		res >> xx >> xy >> xz >> yx >> yy >> yz >> zx >> zy >> zz >> x >> y >> z;
		xyzMatrix< Real > r(xyzMatrix< Real >::rows(
			xx, xy, xz, yx, yy, yz, zx, zy, zz));
		xyzVector< Real > t(x, y, z);
		pose.set_jump(jump_id, Jump(RT(r,t)));
	}

	TR.Debug << "Jumps loaded" << std::endl;
}

void
PoseConformationFeatures::load_chain_endings(
	sessionOP db_session,
	StructureID struct_id,
	Pose & pose
){

	if ( !basic::database::table_exists(db_session, "chain_endings") ) {
		TR.Warning << "chain_endings table does not exist and thus respective data will not be added to the pose!" << std::endl;
		return;
	}

	//note the Conformation object sorts the chain_endings after they are passed in.

	std::string statement_string =
		"SELECT\n"
		"\tend_pos\n"
		"FROM\n"
		"\tchain_endings\n"
		"WHERE\n"
		"\tchain_endings.struct_id=?;";

	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);


	result res(basic::database::safely_read_from_database(stmt));

	vector1< Size > chain_endings;
	while ( res.next() ) {
		Size end_pos;
		res >> end_pos;
		chain_endings.push_back(end_pos);
	}
	pose.conformation().chain_endings(chain_endings);

	TR.Debug << "Chain endings loaded" << std::endl;
}

std::string PoseConformationFeatures::type_name() const {
	return class_name();
}

std::string PoseConformationFeatures::class_name() {
	return "PoseConformationFeatures";
}

void PoseConformationFeatures::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::features::xsd_type_definition_w_attributes( xsd, class_name(), "Records all the chain endings (the final residue of each chain, including the final residue of the Pose overall) as features", attlist );
}

std::string PoseConformationFeaturesCreator::type_name() const {
	return PoseConformationFeatures::class_name();
}

protocols::features::FeaturesReporterOP
PoseConformationFeaturesCreator::create_features_reporter() const {
	return protocols::features::FeaturesReporterOP( new PoseConformationFeatures );
}

void PoseConformationFeaturesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PoseConformationFeatures::provide_xml_schema( xsd );
}



} // namespace
} // namespace
