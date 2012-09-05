// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license.
// (c) The Rosetta software is developed by the contributing members of the
// (c) Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org.
// (c) Questions about this can be addressed to University of Washington UW
// (c) TechTransfer, email: license@u.washington.edu.

/// @file   tester.cc
/// @brief  This is simply a generic pilot app for testing changes.
/// @author Labonte

// includes
#include <iostream>
//#include <algorithm>

#include <basic/Tracer.hh>
#include <devel/init.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
//#include <core/pose/PDBInfo.hh>
//#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
//#include <core/kinematics/MoveMap.hh>

//#include <utility/vector1.hh>
//#include <utility/vector0.hh>

int main(int argc, char *argv[])
{
	using namespace std;
	using namespace core;
	using namespace import_pose;
	using namespace pose;
	//using namespace utility;
	//using namespace kinematics;

	// initialize core
	devel::init(argc, argv);

	// declare variables
	Pose test_pose;

	// import a test pose
	pose_from_pdb(test_pose, "/home/labonte/Workspace/test_input/test.pdb");


	// test new set_chi_true_range() method
	/*MoveMap mm = MoveMap();
	mm.set_chi_true_range(5, 10);

	mm.show(20);*/


	// test new vectorL methods
	/*vector1<char> my_map;
	vector0<char> my_other_map;

	my_map.push_back('A');
	my_map.push_back('B');
	my_map.push_back('C');

	cout << my_map[1] << endl;
	cout << my_map[2] << endl;
	cout << my_map[3] << endl;

	if (my_map.contains('B')) {
		cout << "found!" << endl;
	}
	else {
		cout << "not found!" << endl;
	}
	cout << "B was found at " << my_map.index_of('B') << endl;

	if (my_map.contains('D')) {
		cout << "found!" << endl;
	}
	else {
		cout << "not found!" << endl;
	}
	cout << "D was found at " << my_map.index_of('D') << endl;

	my_other_map.push_back('E');
	my_other_map.push_back('F');
	my_other_map.push_back('G');

	cout << my_other_map[0] << endl;
	cout << my_other_map[1] << endl;
	cout << my_other_map[2] << endl;

	if (my_other_map.contains('F')) {
		cout << "found!" << endl;
	}
	else {
		cout << "not found!" << endl;
	}
	cout << "F was found at " << my_other_map.index_of('F') << endl;

	if (my_other_map.contains('D')) {
		cout << "found!" << endl;
	}
	else {
		cout << "not found!" << endl;
	}
	cout << "D was found at " << my_other_map.index_of('D') << endl;

	cout << "##########" << endl;*/


	// test new update_pose_chains_from_pdb_chains() method
	/*for (Size i = 50; i <= 100; ++i) {
		test_pose.pdb_info()->chain(i, 'X');
	}

	test_pose.update_pose_chains_from_pdb_chains();

	// test modified split_by_chain() method
	Pose new_pose;

	vector1<pose::PoseOP> pose_collection;

	pose_collection = test_pose.split_by_chain();

	new_pose = *pose_collection[2];
	new_pose.dump_pdb("/home/labonte/Code/split.pdb", "");

	test_pose.dump_pdb("/home/labonte/Code/output.pdb", "");*/


	// test atom_name()

	Pose::Residue res2 = test_pose.residue(2);

	cout << res2.atom_name(0) << endl;
	//cout << res2.atom_type(0) << endl;
	cout << res2.atom_type_index(0) << endl;
}




