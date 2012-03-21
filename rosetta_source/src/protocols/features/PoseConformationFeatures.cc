// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/PoseConformationFeatures.cc
/// @brief  report comments stored with each pose
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/PoseConformationFeatures.hh>

//External
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>

// Project Headers
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/io/pdb/file_data.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/types.hh>
#include <utility/sql_database/PrimaryKey.hh>
#include <utility/sql_database/ForeignKey.hh>
#include <utility/sql_database/Column.hh>
#include <utility/sql_database/Schema.hh>
#include <utility/sql_database/Constraint.hh>


// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/database/sql_utils.hh>


// Numeric Headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

// Boost Headers
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <cmath>
#include <sstream>

namespace protocols{
namespace features{

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
using core::chemical::ResidueTypeSetCAP;
using core::chemical::ChemicalManager;
using core::chemical::FA_STANDARD;
using core::chemical::CENTROID;
using core::io::pdb::pose_from_pose;
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

string
PoseConformationFeatures::type_name() const { return "PoseConformationFeatures"; }

string
PoseConformationFeatures::schema() const {
    using namespace utility::sql_database;
    
    Column struct_id("struct_id",DbUUID(), false /*not null*/, false /*don't autoincrement*/);
    
    /******pose_conformations******/
    Schema pose_conformations("pose_conformations");
    pose_conformations.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));

    pose_conformations.add_column( Column("annotated_sequence", DbText()) );
    pose_conformations.add_column( Column("total_residue", DbInteger()) );
    pose_conformations.add_column( Column("fullatom", DbBoolean()) );
        
    
    /******pose_trees******/
    Schema fold_trees("fold_trees");
    fold_trees.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));
    
    fold_trees.add_column( Column("start_res", DbInteger()) );
    fold_trees.add_column( Column("start_atom", DbText()) );
    fold_trees.add_column( Column("stop_res", DbInteger()) );
    fold_trees.add_column( Column("stop_atom", DbText()) );
    fold_trees.add_column( Column("label", DbInteger()) );
    fold_trees.add_column( Column("keep_stub_in_residue", DbBoolean()) );

    
    /******jumps******/
    Schema jumps("jumps");
    jumps.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));
    jumps.add_column( Column("jump_id", DbInteger()) );
    jumps.add_column( Column("xx", DbInteger()) );
    jumps.add_column( Column("xy", DbInteger()) );
    jumps.add_column( Column("xz", DbInteger()) );
    jumps.add_column( Column("yx", DbInteger()) );
    jumps.add_column( Column("yy", DbInteger()) );
    jumps.add_column( Column("yz", DbInteger()) );
    jumps.add_column( Column("zx", DbInteger()) );
    jumps.add_column( Column("zy", DbInteger()) );
    jumps.add_column( Column("zz", DbInteger()) );
    jumps.add_column( Column("x", DbInteger()) );
    jumps.add_column( Column("y", DbInteger()) );
    jumps.add_column( Column("z", DbInteger()) );

    /******chain_endings******/
    Schema chain_endings("chain_endings");
    chain_endings.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));
    chain_endings.add_column( Column("end_pos", DbInteger()) );

    return pose_conformations.print() + "\n" + fold_trees.print() + "\n" + jumps.print() + "\n" + chain_endings.print();
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
	boost::uuids::uuid struct_id,
	sessionOP db_session
){
	vector1< Size > residue_indices;
	for(Size i = 1; i <= relevant_residues.size(); ++i){
		if(relevant_residues[i]) residue_indices.push_back(i);
	}
	Pose* pose; // I had serious memory corruption when I tried to make this an owning pointer
	if (residue_indices.size() == pose_orig.n_residue()){
		pose = const_cast<Pose * >(&pose_orig);
	}
	else{
		pose_from_pose( *pose, pose_orig, residue_indices);
	}

	FoldTree const & fold_tree(pose->conformation().fold_tree());
	//assume non-trivial fold_tree only if more than one edge, i.e., EDGE 1 <nres> -1
	//cppdb::transaction transact_guard(*db_session);
    
	std::string fold_tree_string = "INSERT INTO fold_trees (struct_id, start_res, start_atom, stop_res, stop_atom, label, keep_stub_in_residue) VALUES (?,?,?,?,?,?,?);";
	statement fold_tree_statement(basic::database::safely_prepare_statement(fold_tree_string,db_session));
	for (FoldTree::const_iterator
			it = fold_tree.begin(), it_end = fold_tree.end(); it != it_end; ++it) {
		int start_res(it->start()), stop_res(it->stop()), label(it->label());
		string start_atom(it->start_atom()), stop_atom(it->stop_atom());
		bool keep_stub_in_residue(it->keep_stub_in_residue());

		fold_tree_statement.bind(1,struct_id);
		fold_tree_statement.bind(2,start_res);
		fold_tree_statement.bind(3,start_atom);
		fold_tree_statement.bind(4,stop_res);
		fold_tree_statement.bind(5,stop_atom);
		fold_tree_statement.bind(6,label);
		fold_tree_statement.bind(7,keep_stub_in_residue);
		basic::database::safely_write_to_database(fold_tree_statement);

	}

	std::string jump_string = "INSERT INTO jumps (struct_id, jump_id, xx, xy, xz, yx, yy, yz, zx, zy, zz, x, y, z) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?);";
	statement jump_statement(basic::database::safely_prepare_statement(jump_string,db_session));
 	for (Size nr = 1; nr <= fold_tree.num_jump(); nr++)  {
		Jump const & jump(pose->jump(nr));
		xyzMatrix< Real > const & r(jump.get_rotation());
		Real xx(r.xx()), xy(r.xy()), xz(r.xz());
		Real yx(r.yx()), yy(r.yy()), yz(r.yz());
		Real zx(r.zx()), zy(r.zy()), zz(r.zz());
		Vector const & t(jump.get_translation());
		Real x(t.x()), y(t.y()), z(t.z());

		jump_statement.bind(1,struct_id);
		jump_statement.bind(2,nr);
		jump_statement.bind(3,xx);
		jump_statement.bind(4,xy);
		jump_statement.bind(5,xz);
		jump_statement.bind(6,yx);
		jump_statement.bind(7,yy);
		jump_statement.bind(8,yz);
		jump_statement.bind(9,zx);
		jump_statement.bind(10,zy);
		jump_statement.bind(11,zz);
		jump_statement.bind(12,x);
		jump_statement.bind(13,y);
		jump_statement.bind(14,z);
		basic::database::safely_write_to_database(jump_statement);
	}

 	std::string chain_ending_string = "INSERT INTO chain_endings (struct_id, end_pos) VALUES (?,?);";
 	statement chain_ending_statement(basic::database::safely_prepare_statement(chain_ending_string,db_session));
	foreach(Size end_pos, pose->conformation().chain_endings()){

		chain_ending_statement.bind(1,struct_id);
		chain_ending_statement.bind(2,end_pos);
		basic::database::safely_write_to_database(chain_ending_statement);

	}

	bool ideal = true;
	core::conformation::Conformation const & conformation(pose->conformation());
	for(core::Size resn = 1; resn <= pose->n_residue();++resn)
	{
		bool residue_status(core::conformation::is_ideal_position(resn,conformation));
		if(!residue_status)
		{
			ideal = false;
			break;
		}
	}

	string annotated_sequence(pose->annotated_sequence(true));

	std::string pose_conformation_string = "INSERT INTO pose_conformations (struct_id, annotated_sequence, total_residue, fullatom) VALUES (?,?,?,?);";
	statement pose_conformation_statement(basic::database::safely_prepare_statement(pose_conformation_string,db_session));

	pose_conformation_statement.bind(1,struct_id);
	pose_conformation_statement.bind(2,annotated_sequence);
	pose_conformation_statement.bind(3,pose->total_residue());
	pose_conformation_statement.bind(4,pose->is_fullatom());

	basic::database::safely_write_to_database(pose_conformation_statement);

	//transact_guard.commit();
	return 0;
}

void PoseConformationFeatures::delete_record(
	boost::uuids::uuid struct_id,
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
	boost::uuids::uuid struct_id,
	Pose & pose
){
	load_sequence(db_session, struct_id, pose);
	load_fold_tree(db_session, struct_id, pose);
	load_jumps(db_session, struct_id, pose);
	load_chain_endings(db_session, struct_id, pose);
}

void
PoseConformationFeatures::load_sequence(
	sessionOP db_session,
	boost::uuids::uuid struct_id,
	Pose & pose
){

	std::string statement_string =
		"SELECT\n"
		"	annotated_sequence,\n"
		"	total_residue,\n"
		"	fullatom\n"
		"FROM\n"
		"	pose_conformations\n"
		"WHERE\n"
		"	pose_conformations.struct_id = ?;";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(stmt));

	if(!res.next()){
		stringstream error_message;
		error_message << "Unable to locate structure with struct_id '" << to_string(struct_id) << "'";
		utility_exit_with_message(error_message.str());
	}
	string annotated_sequence;
	Size total_residue, fullatom;
	res >> annotated_sequence >> total_residue >> fullatom;

	ResidueTypeSetCAP residue_set(ChemicalManager::get_instance()->residue_type_set(
		fullatom ? FA_STANDARD : CENTROID));
	make_pose_from_sequence(pose, annotated_sequence, *residue_set);
	runtime_assert(pose.total_residue() == total_residue );

}


void
PoseConformationFeatures::load_fold_tree(
	sessionOP db_session,
	boost::uuids::uuid struct_id,
	Pose & pose
){

	statement stmt = (*db_session) <<
				"SELECT\n"
				"	start_res,\n"
				"	start_atom,\n"
				"	stop_res,\n"
				"	stop_atom,\n"
				"	label,\n"
				"	keep_stub_in_residue\n"
				"FROM\n"
				"	fold_trees\n"
				"WHERE\n"
				"	fold_trees.struct_id=?;" << struct_id;

	result res(basic::database::safely_read_from_database(stmt));

	FoldTree t = FoldTree();
	while(res.next()){
		int start_res, stop_res, label;
		string start_atom, stop_atom;
		int keep_stub_in_residue;
		res >> start_res >> start_atom >> stop_res >> stop_atom >> label >> keep_stub_in_residue;
		if(label == -2 || label > 0){ //CHEMICAL or JUMP
			t.add_edge(Edge(
				start_res, stop_res, label, start_atom, stop_atom, keep_stub_in_residue));
		} else {
			t.add_edge(Edge(start_res, stop_res, label, "", "", keep_stub_in_residue));
		}
	}
	// TODO verify that pose.fold_tree(t) is ok (not cleared from the stack)
	pose.fold_tree(t);
}

void
PoseConformationFeatures::load_jumps(
	sessionOP db_session,
	boost::uuids::uuid struct_id,
	Pose & pose
){
	//note the Conformation object sorts the chain_endings after they are passed in.
	std::string db_mode(basic::options::option[basic::options::OptionKeys::inout::database_mode]);

	statement stmt;
	if(db_mode == "sqlite3")
	{
		std::string statement_string =
			"SELECT\n"
			"	jump_id,\n"
			"	xx REAL,\n"
			"	xy REAL,\n"
			"	xz REAL,\n"
			"	yx REAL,\n"
			"	yy REAL,\n"
			"	yz REAL,\n"
			"	zx REAL,\n"
			"	zy REAL,\n"
			"	zz REAL,\n"
			"	x REAL,\n"
			"	y REAL,\n"
			"	z REAL\n"
			"FROM\n"
			"	jumps\n"
			"WHERE\n"
			"	jumps.struct_id=?;";
		stmt = basic::database::safely_prepare_statement(statement_string,db_session);

	}else if(db_mode == "mysql")
	{

		std::string statement_string =
			"SELECT\n"
			"	jump_id,\n"
			"	xx ,\n"
			"	xy ,\n"
			"	xz ,\n"
			"	yx ,\n"
			"	yy ,\n"
			"	yz ,\n"
			"	zx ,\n"
			"	zy ,\n"
			"	zz ,\n"
			"	x ,\n"
			"	y ,\n"
			"	z \n"
			"FROM\n"
			"	jumps\n"
			"WHERE\n"
			"	jumps.struct_id=?;";
		stmt = basic::database::safely_prepare_statement(statement_string,db_session);

	}else
	{
		utility_exit_with_message("the database mode needs to be 'mysql' or 'sqlite3'");
	}

	stmt.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(stmt));
	while(res.next()){
		Size jump_id;
		Real xx, xy, xz, yx, yy, yz, zx, zy, zz, x, y, z;
		res >> jump_id;
		res >> xx >> xy >> xz >> yx >> yy >> yz >> zx >> zy >> zz >> x >> y >> z;
		xyzMatrix< Real > r(xyzMatrix< Real >::rows(
			xx, xy, xz, yx, yy, yz, zx, zy, zz));
		xyzVector< Real > t(x, y, z);
		pose.set_jump(jump_id, Jump(RT(r,t)));
	}
}

void
PoseConformationFeatures::load_chain_endings(
	sessionOP db_session,
	boost::uuids::uuid struct_id,
	Pose & pose
){

	//note the Conformation object sorts the chain_endings after they are passed in.

	std::string statement_string =
		"SELECT\n"
		"	end_pos\n"
		"FROM\n"
		"	chain_endings\n"
		"WHERE\n"
		"	chain_endings.struct_id=?;";

	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);


	result res(basic::database::safely_read_from_database(stmt));

	vector1< Size > chain_endings;
	while(res.next()){
		Size end_pos;
		res >> end_pos;
		chain_endings.push_back(end_pos);
	}
	pose.conformation().chain_endings(chain_endings);
}



} // namespace
} // namespace
