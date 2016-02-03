// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/pose/copydofs//CopyDofsTest.cxxtest.hh
/// @brief  test CopyDofs by copying internal DOFS from structured pose into new pose, and checking RMSD
/// @author Rhiju Das (rhiju@stanford.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <core/pose/copydofs/CopyDofs.hh>
#include <core/pose/copydofs/util.hh>


// Core Headers
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/MiniPose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/rna/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/import_pose/import_pose.hh>

// Protocol Headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("CopyDofsTest");


class CopyDofsTest : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();

	}

	void tearDown(){

	}



	void test_copy_dofs_into_extended(){

		using namespace core;
		using namespace core::pose;
		using namespace core::pose::rna;
		using namespace core::pose::copydofs;
		using namespace core::import_pose;
		using namespace core::scoring;
		using namespace core::kinematics;

		PoseOP pose_op = pose_from_file( "core/pose/copydofs/excise_bp_7_14_min_again_relax_full_RNA_P_overlap_reR-hbond_sc_wildtype_bound.pdb" , core::import_pose::PDB_file); // mixed RNA, protein, cutpoints.
		Pose & pose = *pose_op;
		figure_out_reasonable_rna_fold_tree( pose );
		virtualize_5prime_phosphates( pose );

		Pose new_pose;
		make_pose_from_sequence( new_pose, pose.sequence(), pose.residue_type(1).residue_type_set(), false /*auto_termin*/ );

		// following fold_tree fix up is necessary unfortunately -- carried out in protocols/stepwise/modeler/util.cc,
		// copying code here to keep this test internal to core.
		FoldTree f( pose.fold_tree() );
		f.set_jump_atoms( 1, " CA ", " CA ", true /*keep stub in residue*/ );
		f.reassign_atoms_for_intra_residue_stubs();
		new_pose.fold_tree( f );
		TS_ASSERT( !f.is_simple_tree() );
		//  new_pose.dump_pdb( "extended.pdb" );

		TS_ASSERT_DIFFERS( pose.annotated_sequence(), new_pose.annotated_sequence() );

		// following 'writes out' what happens in core::pose::copydofs::copy_dofs_match_atom_names()
		std::map < core::id::AtomID , core::id::AtomID > atom_id_map;
		MiniPose const reference_pose( pose );
		std::map< Size, Size > res_map;
		for ( Size n = 1; n <= pose.total_residue(); ++n ) res_map[n] = n;
		setup_atom_id_map_match_atom_names( atom_id_map, res_map, new_pose, reference_pose );

		// copies over all atoms (except virtuals) from reference_pose to extended pose.
		CopyDofs copy_dofs( reference_pose, atom_id_map );
		copy_dofs.apply( new_pose );

		TS_ASSERT( rms_at_corresponding_atoms_no_super( new_pose, pose, atom_id_map ) > 10.0 );

		// really non-trivial test -- PDBs must be superimposable now.
		TS_ASSERT_LESS_THAN( superimpose_pose( new_pose, pose, atom_id_map ), 1.0e-4 );
		// above did the superposition.
		TS_ASSERT_LESS_THAN( rms_at_corresponding_atoms_no_super( new_pose, pose, atom_id_map ), 1.0e-4 );

		//  pose.dump_pdb( "test.pdb" );
		//  new_pose.dump_pdb( "copydofs_superimpose.pdb" );
	}



};



