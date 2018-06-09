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
#include <protocols/antibody/snugdock/SnugDockProtocol.hh>
#include <protocols/antibody/AntibodyInfo.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.antibody.snugdock.SnugDockProtocol_setup_ab_ag_foldtree");
class SnugDockProtocol_setup_ab_ag_foldtree_tests : public CxxTest::TestSuite
{
public:
	void setUp() {
		// shared initialization goes here
		// do we set -s here? (example below)
		// core_init_with_additional_options( "-antibody:exclude_homologs true -antibody:exclude_homologs_fr_cutoff 60 -out:level 500 -antibody:n_multi_templates 3 -antibody:ocd_cutoff 5" );
		core_init_with_additional_options( "-input_ab_scheme AHO_Scheme" );
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

	void test_setup_ab_ag_foldtree() {
		// load pose and setup foldtree -- first pose is just VH-VL-Ag, second is VH-Ag_C-Ag_D
		core::pose::Pose antibody_with_antigen;
		core::import_pose::pose_from_file(antibody_with_antigen, "protocols/antibody/aho_with_antigen.pdb", core::import_pose::PDB_file);
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
		// Ab_COM is 6.784177419354839 67.59274193548387 65.39931451612900
		compareVectors( antibody_with_antigen.xyz( core::id::AtomID(1, 1) ),
			numeric::xyzVector<core::Real>(6.78417, 67.59274, 65.39931),
			0.0001 );

		// Ag_COM is 28.96556470588236 73.00931764705884 83.78094117647058
		compareVectors( antibody_with_antigen.xyz( core::id::AtomID(1, 2) ),
			numeric::xyzVector<core::Real>(28.96556, 73.00931, 83.78094),
			0.0001 );

		// H_COM is 8.106810218978101 74.34232116788324 59.3131678832116
		compareVectors( antibody_with_antigen.xyz( core::id::AtomID(1, 3) ),
			numeric::xyzVector<core::Real>(8.10681, 74.34232, 59.31316),
			0.0001 );

		// L_COM is 5.151738738738739 59.26218018018018 72.91104504504503
		compareVectors( antibody_with_antigen.xyz( core::id::AtomID(1, 4) ),
			numeric::xyzVector<core::Real>(5.15173, 59.26218, 72.911104),
			0.0001 );

		// S_COM is 28.96556470588236 73.00931764705884 83.78094117647058
		compareVectors( antibody_with_antigen.xyz( core::id::AtomID(1, 5) ),
			numeric::xyzVector<core::Real>(28.96556, 73.00931, 83.78094),
			0.0001 );

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
		// Ab_COM is 37.89114960629922 8.466818897637795 5.953889763779528
		compareVectors( antibody_with_antigen.xyz( core::id::AtomID(1, 1) ),
			numeric::xyzVector<core::Real>(37.89114, 8.46681, 5.95388),
			0.0001 );

		// Ag_COM is 10.51437500000000 9.730505952380947 8.777624999999995
		compareVectors( antibody_with_antigen.xyz( core::id::AtomID(1, 2) ),
			numeric::xyzVector<core::Real>(10.51437, 9.73050, 8.77762),
			0.0001 );

		// H_COM is 37.89114960629922 8.466818897637795 5.953889763779528
		compareVectors( antibody_with_antigen.xyz( core::id::AtomID(1, 3) ),
			numeric::xyzVector<core::Real>(37.89114, 8.46681, 5.95388),
			0.0001 );

		// C_COM is 6.765279069767443 15.94729069767442 20.96729069767441
		compareVectors( antibody_with_antigen.xyz( core::id::AtomID(1, 4) ),
			numeric::xyzVector<core::Real>(6.76527, 15.94729, 20.96729),
			0.0001 );

		// D_COM is 14.44635365853658 3.210463414634147 -4.006658536585367
		compareVectors( antibody_with_antigen.xyz( core::id::AtomID(1, 5) ),
			numeric::xyzVector<core::Real>(14.44635, 3.21046, -4.00665),
			0.0001 );

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
