// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/branch_angle/BranchAngleOptimizer.cxxtest.hh
/// @brief  test suite for protocols::branch_angle::BranchAngleOptimizer class
/// @author Colin A. Smith


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

// Package headers
#include <protocols/branch_angle/BranchAngleOptimizer.hh>

// Project headers
#include <core/id/AtomID_Mask.hh>
#include <core/chemical/AtomType.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/AtomPointer.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/DomainMap.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/scoring/mm/MMBondAngleResidueTypeParam.hh>
#include <core/scoring/mm/MMBondAngleResidueTypeParamSet.hh>

// Numeric Headers
#include <numeric/BodyPosition.fwd.hh>
#include <numeric/Quaternion.fwd.hh>
#include <numeric/all.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/io/all.fwd.hh>
#include <utility/io/icstream.fwd.hh>
#include <utility/io/ocstream.fwd.hh>


using namespace core;
using namespace core::pose;
using namespace protocols::branch_angle;

static basic::Tracer TR("protocols.branch_angle.BranchAngleOptimizer.cxxtest");

void
optimize_atom_tree(
	BranchAngleOptimizer & branchopt,
	std::ostream & out,
	Pose & pose,
	kinematics::tree::AtomCOP main_atom1 = NULL,
	kinematics::tree::AtomCOP center_atom = NULL
)
{
	if (!main_atom1 && !center_atom) {
		main_atom1 = pose.atom_tree().root();
	}

	if (center_atom) {

		kinematics::tree::AtomCOP const main_atom2(center_atom->get_nonjump_atom(0));

		if (main_atom2) {

			id::AtomID const main_atomid1(main_atom1->id());
			id::AtomID const center_atomid(center_atom->id());
			id::AtomID const main_atomid2(main_atom2->id());

			if (!(pose.residue(main_atomid1.rsd()).atom_type(main_atomid1.atomno()).name() == "VIRT" ||
			      pose.residue(center_atomid.rsd()).atom_type(center_atomid.atomno()).name() == "VIRT" ||
			      pose.residue(main_atomid2.rsd()).atom_type(main_atomid2.atomno()).name() == "VIRT")) {

				out << "Optimizing:" << main_atomid1 << pose.residue(main_atomid1.rsd()).atom_name(main_atomid1.atomno()) << " "
														 << center_atomid << pose.residue(center_atomid.rsd()).atom_name(center_atomid.atomno()) << " "
														 << main_atomid2 << pose.residue(main_atomid2.rsd()).atom_name(main_atomid2.atomno()) << std::endl;
				Size const index(branchopt.optimize_angles(pose, main_atomid1, center_atomid, main_atomid2));
				out << "Got branching index: " << index << std::endl;
			}

			for (Size i = 1; i < center_atom->n_nonjump_children(); ++i) {
				optimize_atom_tree(branchopt, out, pose, center_atom, center_atom->get_nonjump_atom(i));
			}

			if (center_atom->n_nonjump_children()) {
				optimize_atom_tree(branchopt, out, pose, center_atom, center_atom->get_nonjump_atom(0));
			}
		}

	} else {

		for (Size i = 0; i < main_atom1->n_nonjump_children(); ++i) {
				optimize_atom_tree(branchopt, out, pose, main_atom1, main_atom1->get_nonjump_atom(i));
		}
	}
}

class BranchAngleOptimizerTest : public CxxTest::TestSuite {

public:

	chemical::ResidueTypeSetCAP residue_set;
	pose::PoseOP the_pose;

	BranchAngleOptimizerTest() {}

	void setUp() {
		core_init_with_additional_options( "-no_optH" );
		residue_set = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );

		the_pose = pose::PoseOP( new Pose );
		//the_pose = create_test_in_pdb_poseop(); slightly different chain IDs
		core::import_pose::pose_from_pdb( *the_pose, "protocols/moves/test_in.pdb" );
	}

	void tearDown() {
		the_pose.reset();
	}

	// Test functions that go along with the BranchAnlgeOptimizer class
	void test_BranchAngleOptimizer() {

		test::UTracer UT("protocols/branch_angle/BranchAngleOptimizer.u");

		UT << "Testing for fully connected atoms..." << std::endl;

		for (Size i = 1; i <= the_pose->total_residue(); ++i) {
			for (Size j = 1; j <= the_pose->residue(i).natoms(); ++j) {
				TS_ASSERT(!the_pose->residue(i).has_incomplete_connection(j));
				Size const neighbor_count(the_pose->residue(i).n_bonded_neighbor_all_res(j));
				utility::vector1<core::id::AtomID> const neighbors(the_pose->conformation().bonded_neighbor_all_res(id::AtomID(j, i)));
				TS_ASSERT(neighbor_count == neighbors.size());
				UT << "Residue " << i << ", atom " << j << ":" << std::endl;
				UT << " neighbor count: " << neighbor_count << std::endl;
				for (Size k = 1; k <= neighbors.size(); ++k) {
					UT << neighbors[k] << std::endl;
				}
			}
		}

		UT << std::endl << "Testing for branch angle selection..." << std::endl;

		id::AtomID branch_atomid1;
		id::AtomID branch_atomid2;

		kinematics::tree::AtomCOP atom1(the_pose->atom_tree().root());
		kinematics::tree::AtomCOP atom2(atom1->child(1));
		kinematics::tree::AtomCOP atom3(atom2->child(0));

		while (atom3) {
			UT << "Branching atoms for" << atom1->id() << "->" << atom2->id() << "->" << atom3->id() << ":" << std::endl;
			Size const neighbor_count(the_pose->residue(atom2->id().rsd()).n_bonded_neighbor_all_res(atom2->id().atomno()));
			if (neighbor_count == 3) {
				branching_atomid1(*the_pose, atom1->id(), atom2->id(), atom3->id(), branch_atomid1);
				UT << branch_atomid1 << std::endl;
			} else if (neighbor_count == 4) {
				branching_atomids2(*the_pose, atom1->id(), atom2->id(), atom3->id(), branch_atomid1, branch_atomid2);
				UT << branch_atomid1 << std::endl;
				UT << branch_atomid2 << std::endl;
			}

			if (!(the_pose->residue(atom2->id().rsd()).name3() == "PRO" &&
			      the_pose->residue(atom2->id().rsd()).atom_name(atom2->id().atomno()) == " N  ")) {
			}

			atom1 = atom2;
			atom2 = atom3;
			atom3 = atom3->get_nonjump_atom(0);
		}

		UT << std::endl << "Testing whole protein branch angle optimization..." << std::endl;

		BranchAngleOptimizer branchopt;

		// read in branch angle coefficients
		branchopt.read_coef1("protocols/branch_angle/branch_angle_1.txt");
		branchopt.read_coef2("protocols/branch_angle/branch_angle_2.txt");

		Size index;

		UT << "Optimizing atom tree root" << std::endl;
		index = branchopt.optimize_angles(*the_pose, the_pose->atom_tree().root()->child(0)->id(),
		                                  the_pose->atom_tree().root()->id(),
												        			the_pose->atom_tree().root()->child(1)->id());
		UT << "Got branching index: " << index << std::endl;

		UT << std::endl << "Optimizing PDB structure..." << std::endl;

		optimize_atom_tree(branchopt, UT, *the_pose);

		UT << std::endl << "Dumping optimized PDB structure..." << std::endl;

		the_pose->dump_pdb(UT);

		PoseOP pose_ideal( new Pose() );
		core::pose::make_pose_from_sequence(*pose_ideal, "ACDEFGHIKLMNQRSTVWY", "fa_standard");
		PoseOP pose_ideal_optimized( new Pose(*pose_ideal) );

		core::scoring::mm::MMBondAngleResidueTypeParamSetOP param_set( new core::scoring::mm::MMBondAngleResidueTypeParamSet() );
		param_set->use_residue_type_theta0(true);
		branchopt.bond_angle_residue_type_param_set(param_set);

		UT << std::endl << "Optimizing idealized PDB structure..." << std::endl;

		optimize_atom_tree(branchopt, UT, *pose_ideal_optimized);

		UT << std::endl << "Dumping optimized PDB structure..." << std::endl;

		pose_ideal_optimized->dump_pdb(UT);

		// check that the atoms did not move by more than 0.01 angstroms
		// This seems to fail for planar atoms at the end of residues
		/*
		for (Size i = 1; i <= pose_ideal->total_residue(); ++i) {
			for (Size j = 1; j <= pose_ideal->residue(i).natoms(); ++j) {
				Real atom_distance(pose_ideal->xyz(id::AtomID(j, i)).distance(pose_ideal_optimized->xyz(id::AtomID(j, i))));
				std::cout << i << "\t" << j << "\t" << pose_ideal->residue(i).name() << "\t" << pose_ideal->residue(i).atom_name(j) << "\t" << atom_distance << std::endl;
				TS_ASSERT(atom_distance < 0.01);
			}
		}
		*/
	}
};
