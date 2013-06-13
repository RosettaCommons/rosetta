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
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/MoveMap.hh>

#include <utility/excn/Exceptions.hh>
//#include <utility/vector1.hh>
//#include <utility/vector0.hh>

#include <protocols/simple_moves/BackboneMover.hh>

int main(int argc, char *argv[])
{
    try {

	using namespace std;
	using namespace core;
	using namespace import_pose;
	using namespace pose;
	using namespace utility;
	using namespace kinematics;
	using namespace protocols;

	// initialize core
	devel::init(argc, argv);

	// declare variables
	Pose pose;

	// import a test pose
	//pose_from_pdb(pose, "/home/labonte/Workspace/test_input/test.pdb");

	pose_from_pdb(pose, "/home/labonte/Workspace/Carbohydrates/maltotriose.pdb");

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


	// test carbohydrate methods
	/*cout << "Start: Psi of residue 2: " << pose.psi(2) << endl;
	pose.set_psi(2, 75.0);
	cout << "End: Psi of residue 2: " << pose.psi(2) << endl;*/

	/*cout << "Start: Phi of residue 2: " << pose.phi(2) << endl;
	pose.set_phi(2, 95.0);
	cout << "End: Phi of residue 2: " << pose.phi(2) << endl;*/

	Size n_res = pose.total_residue();

	MoveMapOP mm = new MoveMap();
	mm->set_bb_true_range(2, n_res);

	simple_moves::SmallMover mover;

	mover.movemap(mm);

	cout << "BEFORE MOVE" << endl;
	cout << "Residue: 1" << endl;

	cout << " Chis:" << endl;
	cout << "  1 " << pose.chi(1, 1) << endl;
	cout << "  2 " << pose.chi(2, 1) << endl;
	cout << "  3 " << pose.chi(3, 1) << endl;
	cout << "  4 " << pose.chi(4, 1);
	cout << " (should equal psi(n+1): " << pose.psi(2) << ")" << endl;
	cout << "  5 " << pose.chi(5, 1) << endl;

	for (Size i = 2; i <= n_res - 1; ++i) {
		cout << "Residue: " << i << endl;
		cout << " Phi: " << pose.phi(i) << endl;
		cout << " Psi: " << pose.psi(i) << endl;

		cout << " Chis:" << endl;
		cout << "  1 " << pose.chi(1, i);
		cout << " (should equal phi: " << pose.phi(i) << ")" << endl;
		cout << "  2 " << pose.chi(2, i) << endl;
		cout << "  3 " << pose.chi(3, i) << endl;
		cout << "  4 " << pose.chi(4, i);
		cout << " (should equal psi(n+1): " << pose.psi(i + 1) << ")" << endl;
		cout << "  5 " << pose.chi(5, i) << endl;
	}

	cout << "Residue: " << n_res << endl;
	cout << " Phi: " << pose.phi(n_res) << endl;
	cout << " Psi: " << pose.psi(n_res) << endl;

	cout << " Chis:" << endl;
	cout << "  1 " << pose.chi(1, n_res);
	cout << " (should equal phi: " << pose.phi(n_res) << ")" << endl;
	cout << "  2 " << pose.chi(2, n_res) << endl;
	cout << "  3 " << pose.chi(3, n_res) << endl;
	cout << "  4 " << pose.chi(4, n_res) << endl;
	cout << "  5 " << pose.chi(5, n_res) << endl;

	mover.apply(pose);

	cout << "AFTER MOVE" << endl;

	cout << "Residue: 1" << endl;

	cout << " Chis:" << endl;
	cout << "  1 " << pose.chi(1, 1) << endl;
	cout << "  2 " << pose.chi(2, 1) << endl;
	cout << "  3 " << pose.chi(3, 1) << endl;
	cout << "  4 " << pose.chi(4, 1);
	cout << " (should equal psi(n+1): " << pose.psi(2) << ")" << endl;
	cout << "  5 " << pose.chi(5, 1) << endl;

	for (Size i = 2; i <= n_res - 1; ++i) {
		cout << "Residue: " << i << endl;
		cout << " Phi: " << pose.phi(i) << endl;
		cout << " Psi: " << pose.psi(i) << endl;

		cout << " Chis:" << endl;
		cout << "  1 " << pose.chi(1, i);
		cout << " (should equal phi: " << pose.phi(i) << ")" << endl;
		cout << "  2 " << pose.chi(2, i) << endl;
		cout << "  3 " << pose.chi(3, i) << endl;
		cout << "  4 " << pose.chi(4, i);
		cout << " (should equal psi(n+1): " << pose.psi(i + 1) << ")" << endl;
		cout << "  5 " << pose.chi(5, i) << endl;
	}

	cout << "Residue: " << n_res << endl;
	cout << " Phi: " << pose.phi(n_res) << endl;
	cout << " Psi: " << pose.psi(n_res) << endl;

	cout << " Chis:" << endl;
	cout << "  1 " << pose.chi(1, n_res);
	cout << " (should equal phi: " << pose.phi(n_res) << ")" << endl;
	cout << "  2 " << pose.chi(2, n_res) << endl;
	cout << "  3 " << pose.chi(3, n_res) << endl;
	cout << "  4 " << pose.chi(4, n_res) << endl;
	cout << "  5 " << pose.chi(5, n_res) << endl;

	pose.dump_pdb("/home/labonte/Workspace/test_output/modified_sugar.pdb", "");

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
    }
    return 0;
}
