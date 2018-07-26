// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/antibody/snugdock/SnugDockProtocol_setup_ab_ag_foldtree.cxxtest.hh
/// @brief  test suite for universal FT setup function in SnugDock
/// @author Jeliazko Jeliazkov

#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>

#include <numeric/xyzVector.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <protocols/antibody/snugdock/SnugDockProtocol.hh>
#include <protocols/antibody/AntibodyInfo.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.antibody.snugdock.SnugDockProtocol_setup_ab_ag_foldtree");
class SnugDockProtocol_setup_ab_ag_foldtree_tests : public CxxTest::TestSuite
{
public:
	core::pose::Pose antibody_with_antigen;
	void setUp() {
		// shared initialization goes here
		// do we set -s here? (example below)
		// core_init_with_additional_options( "-antibody:exclude_homologs true -antibody:exclude_homologs_fr_cutoff 60 -out:level 500 -antibody:n_multi_templates 3 -antibody:ocd_cutoff 5" );
		core_init_with_additional_options( "-input_ab_scheme AHO_Scheme" );
		core::import_pose::pose_from_file(antibody_with_antigen, "protocols/antibody/aho_with_antigen.pdb", core::import_pose::PDB_file);
	} //setUp

	void tearDown() {
		// shared finalization goes here
	} //tearDown

	void testEdgeConnectivity(core::kinematics::Edge edge, core::Size start, core::Size stop, core::Size label) {
		TS_ASSERT_EQUALS(edge.start(), start);
		TS_ASSERT_EQUALS(edge.stop(), stop);
		TS_ASSERT_EQUALS(edge.label(), label);
	} //testEdgeConnectivity

	void compareVectors( numeric::xyzVector< core::Real > abc, numeric::xyzVector< core::Real > xyz, core::Real delta) {
		TS_ASSERT_DELTA(abc.x(), xyz.x(), delta);
		TS_ASSERT_DELTA(abc.y(), xyz.y(), delta);
		TS_ASSERT_DELTA(abc.z(), xyz.z(), delta);
	} //compareVectors

	void test_place_VRT_at_residue_COM() {
		core::Size com_resnum = core::pose::residue_center_of_mass( antibody_with_antigen, 1, antibody_with_antigen.size() );
		TR << "COM residue is: " << com_resnum << std::endl;
		// produce VRT with function and compare coordinates
		protocols::antibody::snugdock::SnugDockProtocolOP snugdock( new protocols::antibody::snugdock::SnugDockProtocol() );
		core::conformation::ResidueOP vrt_res = snugdock->place_VRT_at_residue_COM( antibody_with_antigen, 1, antibody_with_antigen.size() );
		// compare relevant coordinates
		compareVectors( antibody_with_antigen.residue(com_resnum).xyz("CA"), vrt_res->xyz("ORIG"), 0.001 );
		compareVectors( antibody_with_antigen.residue(com_resnum).xyz("N"), vrt_res->xyz("X"), 0.001 );
		compareVectors( antibody_with_antigen.residue(com_resnum-1).xyz("C"), vrt_res->xyz("Y"), 0.001 );
	} //test_place_VRT_at_residue_COM()

	void test_setup_ab_ag_foldtree() {
		using namespace core;
		using namespace numeric;
		// load pose and setup foldtree -- first pose is just VH-VL-Ag, second is VH-Ag_C-Ag_D
		protocols::antibody::AntibodyInfoOP ab_info( new protocols::antibody::AntibodyInfo(antibody_with_antigen) );
		protocols::antibody::snugdock::SnugDockProtocolOP snugdock( new protocols::antibody::snugdock::SnugDockProtocol() );
		snugdock->setup_ab_ag_foldtree(antibody_with_antigen, ab_info);

		TR << antibody_with_antigen.fold_tree() << std::endl;

		// test pose for VRTs being prepended: (1) Ab COM, (2) Ag COM, (3) VH COM, (4) VL COM, (5) Ag-A COM
		for ( core::Size i = 1; i <= 5; ++i ) {
			TS_ASSERT(antibody_with_antigen.residue(i).is_virtual_residue());
		}
		// test pose for first protein residue at 6
		TS_ASSERT(!antibody_with_antigen.residue(6).is_virtual_residue());

		// test COM position versus precalculated
		// Ab_COM ORIG is 4.887 71.464 66.427
		// Ab_COM X is    5.550 71.370 67.721
		// Ab_COM Y is    5.430 72.308 68.692
		compareVectors( antibody_with_antigen.residue(1).xyz("ORIG"),
			xyzVector<Real>(4.887, 71.464, 66.427), 0.001);
		compareVectors( antibody_with_antigen.residue(1).xyz("X"),
			xyzVector<Real>(5.550, 71.370, 67.721), 0.001);
		compareVectors( antibody_with_antigen.residue(1).xyz("Y"),
			xyzVector<Real>(5.430, 72.308, 68.692), 0.001);

		// Ag_COM ORIG is 31.365 73.436 85.610
		// Ag_COM X is    30.535 74.606 85.384
		// Ag_COM Y is    30.393 75.130 84.174
		compareVectors( antibody_with_antigen.residue(2).xyz("ORIG"),
			xyzVector<Real>(31.365, 73.436, 85.610), 0.001);
		compareVectors( antibody_with_antigen.residue(2).xyz("X"),
			xyzVector<Real>(30.535, 74.606, 85.384), 0.001);
		compareVectors( antibody_with_antigen.residue(2).xyz("Y"),
			xyzVector<Real>(30.393, 75.130, 84.174), 0.001);


		// H_COM ORIG is 8.869 74.033 59.787
		// H_COM X is    8.022 75.168 60.170
		// H_COM Y is    8.448 76.419 60.318
		compareVectors( antibody_with_antigen.residue(3).xyz("ORIG"),
			xyzVector<Real>(8.869, 74.033, 59.787), 0.001);
		compareVectors( antibody_with_antigen.residue(3).xyz("X"),
			xyzVector<Real>(8.022, 75.168, 60.170), 0.001);
		compareVectors( antibody_with_antigen.residue(3).xyz("Y"),
			xyzVector<Real>(8.448, 76.419, 60.318), 0.001);

		// L_COM ORIG is 4.322 62.794 72.823
		// L_COM X is    4.969 63.972 73.377
		// L_COM Y is    6.288 64.144 73.222
		compareVectors( antibody_with_antigen.residue(4).xyz("ORIG"),
			xyzVector<Real>(4.322, 62.794, 72.823), 0.001);
		compareVectors( antibody_with_antigen.residue(4).xyz("X"),
			xyzVector<Real>(4.969, 63.972, 73.377), 0.001);
		compareVectors( antibody_with_antigen.residue(4).xyz("Y"),
			xyzVector<Real>(6.288, 64.144, 73.222), 0.001);

		// S_COM ORIG is 31.365 73.436 85.610
		// S_COM X is    30.535 74.606 85.384
		// S_COM Y is    30.393 75.130 84.174
		compareVectors( antibody_with_antigen.residue(5).xyz("ORIG"),
			xyzVector<Real>(31.365, 73.436, 85.610), 0.001);
		compareVectors( antibody_with_antigen.residue(5).xyz("X"),
			xyzVector<Real>(30.535, 74.606, 85.384), 0.001);
		compareVectors( antibody_with_antigen.residue(5).xyz("Y"),
			xyzVector<Real>(30.393, 75.130, 84.174), 0.001);

		// test for correct jump connectivity
		utility::vector1< core::kinematics::Edge > jump_edges = antibody_with_antigen.fold_tree().get_jump_edges();
		// edge #1 is a jump 1-2 (Ab-Ag)
		testEdgeConnectivity(jump_edges[1], 1, 2, 1);
		// edge #2 is a jump 2-3 (Ab-H)
		testEdgeConnectivity(jump_edges[2], 1, 3, 2);
		// edge #3 is a jump 3-4 (H-L)
		testEdgeConnectivity(jump_edges[3], 3, 4, 3);
		// edge #4 is a jump 2-5 (Ag-S) -- S is the antigen chain here
		testEdgeConnectivity(jump_edges[4], 2, 5, 4);
		// edge #5 is a jump 3-6 (H-Nter) -- Nter is the start of the corresponding chain
		testEdgeConnectivity(jump_edges[5], 3, 6, 5);
		// edge #6 is a jump 4-143 (L-Nter) [chain H is 137 residues in size]
		testEdgeConnectivity(jump_edges[6], 4, 143, 6);
		// edge #7 is a jump 5-254 (S-Nter) [chain L is 111 residues in size]
		testEdgeConnectivity(jump_edges[7], 5, 254, 7);

		// test CDR loops in FT start is old loop + 5 (n_vrts) - 3 (as in RefineOneCDRLoop)
		// test CDR loops in FT stop is old loop + 5 (n_vrts) + 3 (as in RefineOneCDRLoop)
		// cut point is the old cutpoint + 5
		testEdgeConnectivity( jump_edges[8],
			(ab_info->get_CDR_loop(protocols::antibody::h1)).start() + 5 - 3,
			(ab_info->get_CDR_loop(protocols::antibody::h1)).stop() + 5 + 3,
			8 );
		// H2
		testEdgeConnectivity( jump_edges[9],
			(ab_info->get_CDR_loop(protocols::antibody::h2)).start() + 5 - 3,
			(ab_info->get_CDR_loop(protocols::antibody::h2)).stop() + 5 + 3,
			9 );
		// H3
		testEdgeConnectivity( jump_edges[10],
			(ab_info->get_CDR_loop(protocols::antibody::h3)).start() + 5 - 3,
			(ab_info->get_CDR_loop(protocols::antibody::h3)).stop() + 5 + 3,
			10 );
		// L1
		testEdgeConnectivity( jump_edges[11],
			(ab_info->get_CDR_loop(protocols::antibody::l1)).start() + 5 - 3,
			(ab_info->get_CDR_loop(protocols::antibody::l1)).stop() + 5 + 3,
			11 );
		// L2
		testEdgeConnectivity( jump_edges[12],
			(ab_info->get_CDR_loop(protocols::antibody::l2)).start() + 5 - 3,
			(ab_info->get_CDR_loop(protocols::antibody::l2)).stop() + 5 + 3,
			12 );
		// L3
		testEdgeConnectivity( jump_edges[13],
			(ab_info->get_CDR_loop(protocols::antibody::l3)).start() + 5 - 3,
			(ab_info->get_CDR_loop(protocols::antibody::l3)).stop() + 5 + 3,
			13 );


		// on to the next ab
		core::import_pose::pose_from_file(antibody_with_antigen, "protocols/antibody/4dka.aho.pdb", core::import_pose::PDB_file);
		protocols::antibody::AntibodyInfoOP ab_info_2 ( new protocols::antibody::AntibodyInfo(antibody_with_antigen) );
		protocols::antibody::snugdock::SnugDockProtocolOP snugdock_2( new protocols::antibody::snugdock::SnugDockProtocol() );
		snugdock_2->setup_ab_ag_foldtree(antibody_with_antigen, ab_info_2);

		TR << antibody_with_antigen.fold_tree() << std::endl;

		// test pose for VRTs being prepended: (1) Ab COM, (2) Ag COM, (3) VH COM, (4) Ag-C, (5) Ag-D COM
		for ( core::Size i = 1; i <= 5; ++i ) {
			TS_ASSERT(antibody_with_antigen.residue(i).is_virtual_residue());
		}
		// test pose for first protein residue at 6
		TS_ASSERT(!antibody_with_antigen.residue(6).is_virtual_residue());

		// test COM position versus precalculated
		// Ab_COM ORIG is 36.943 6.785 6.589
		// Ab_COM X is    35.696 7.099 5.896
		// Ab_COM Y is    35.504 6.946 4.587
		compareVectors( antibody_with_antigen.residue(1).xyz("ORIG"),
			xyzVector<Real>(36.943, 6.785, 6.589), 0.001);
		compareVectors( antibody_with_antigen.residue(1).xyz("X"),
			xyzVector<Real>(35.696, 7.099, 5.896), 0.001);
		compareVectors( antibody_with_antigen.residue(1).xyz("Y"),
			xyzVector<Real>(35.504, 6.946, 4.587), 0.001);

		// Ag_COM ORIG is 14.920 12.734 6.892
		// Ag_COM X is    14.650 13.842 7.799
		// Ag_COM Y is    13.408 14.238 8.074
		compareVectors( antibody_with_antigen.residue(2).xyz("ORIG"),
			xyzVector<Real>(14.920, 12.734, 6.892), 0.001);
		compareVectors( antibody_with_antigen.residue(2).xyz("X"),
			xyzVector<Real>(14.650, 13.842, 7.799), 0.001);
		compareVectors( antibody_with_antigen.residue(2).xyz("Y"),
			xyzVector<Real>(13.408, 14.238, 8.074), 0.001);

		// H_COM ORIG is 36.943 6.785 6.589
		// H_COM X is    35.696 7.099 5.896
		// H_COM Y is    35.504 6.946 4.587
		compareVectors( antibody_with_antigen.residue(3).xyz("ORIG"),
			xyzVector<Real>(36.943, 6.785, 6.589), 0.001);
		compareVectors( antibody_with_antigen.residue(3).xyz("X"),
			xyzVector<Real>(35.696, 7.099, 5.896), 0.001);
		compareVectors( antibody_with_antigen.residue(3).xyz("Y"),
			xyzVector<Real>(35.504, 6.946, 4.587), 0.001);

		// C_COM ORIG is 5.678 15.698 19.833
		// C_COM X is    6.550 16.583 20.603
		// C_COM Y is    7.484 16.137 21.436
		compareVectors( antibody_with_antigen.residue(4).xyz("ORIG"),
			xyzVector<Real>(5.678, 15.698, 19.833), 0.001);
		compareVectors( antibody_with_antigen.residue(4).xyz("X"),
			xyzVector<Real>(6.550, 16.583, 20.603), 0.001);
		compareVectors( antibody_with_antigen.residue(4).xyz("Y"),
			xyzVector<Real>(7.484, 16.137, 21.436), 0.001);

		// D_COM ORIG is 13.676 2.121 -3.343
		// D_COM X is    14.803 2.287 -4.257
		// D_COM Y is    14.827 3.217 -5.200
		compareVectors( antibody_with_antigen.residue(5).xyz("ORIG"),
			xyzVector<Real>(13.676, 2.121, -3.343), 0.001);
		compareVectors( antibody_with_antigen.residue(5).xyz("X"),
			xyzVector<Real>(14.803, 2.287, -4.257), 0.001);
		compareVectors( antibody_with_antigen.residue(5).xyz("Y"),
			xyzVector<Real>(14.827, 3.217, -5.200), 0.001);

		// test for correct jump connectivity
		jump_edges = antibody_with_antigen.fold_tree().get_jump_edges();
		// edge #1 is a jump 1-2 (Ab-Ag)
		testEdgeConnectivity(jump_edges[1], 1, 2, 1);
		// edge #2 is a jump 2-3 (Ab-H)
		testEdgeConnectivity(jump_edges[2], 1, 3, 2);
		// edge #3 is a jump 3-4 (Ag-C) -- C is the antigen chain here
		testEdgeConnectivity(jump_edges[3], 2, 4, 3);
		// edge #4 is a jump 2-5 (Ag-D)
		testEdgeConnectivity(jump_edges[4], 4, 5, 4);
		// edge #5 is a jump 3-6 (H-Nter) -- Nter is the start of the corresponding chain
		testEdgeConnectivity(jump_edges[5], 3, 6, 5);
		// edge #6 is a jump 4-133 (C-Nter) [chain H is 127 residues in size]
		testEdgeConnectivity(jump_edges[6], 4, 133, 6);
		// edge #7 is a jump 5-219 (D-Nter) [chain C is 86 residues in size]
		testEdgeConnectivity(jump_edges[7], 5, 219, 7);

		// test CDR loops in FT start is old loop + 5 (n_vrts) - 3 (as in RefineOneCDRLoop)
		// test CDR loops in FT stop is old loop + 5 (n_vrts) + 3 (as in RefineOneCDRLoop)
		// cut point is the old cutpoint + 5
		testEdgeConnectivity( jump_edges[8],
			(ab_info_2->get_CDR_loop(protocols::antibody::h1)).start() + 5 - 3,
			(ab_info_2->get_CDR_loop(protocols::antibody::h1)).stop() + 5 + 3,
			8 );
		// H2
		testEdgeConnectivity( jump_edges[9],
			(ab_info_2->get_CDR_loop(protocols::antibody::h2)).start() + 5 - 3,
			(ab_info_2->get_CDR_loop(protocols::antibody::h2)).stop() + 5 + 3,
			9 );
		// H3
		testEdgeConnectivity( jump_edges[10],
			(ab_info_2->get_CDR_loop(protocols::antibody::h3)).start() + 5 - 3,
			(ab_info_2->get_CDR_loop(protocols::antibody::h3)).stop() + 5 + 3,
			10 );

	}
};
