// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/ResfileReader.cxxtest.hh
/// @brief  test suite for resfile reader
/// @author Steven Lewis

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Project Headers
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/rotamer_trials.hh>

//#include <protocols/simple_moves/PackRotamersMover.hh> //I wish

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh> //CENTROID, FA_STANDARD

// Numeric Headers

// Utility Headers
// AUTO-REMOVED #include <core/init/init.hh>
// AUTO-REMOVED #include <basic/Tracer.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>



// --------------- Test Class --------------- //

class FixHisTautomerTest : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	void setUp() {
		core_init_with_additional_options( "-restore_pre_talaris_2013_behavior -override_rsd_type_limit" );
	}

	void tearDown() {
	}

	// ------------- Helper Functions ------------- //
	///@brief check what happens to a packer task when the only instruction is that it should be repacked
	core::pack::task::PackerTaskOP
	task_repack( core::pose::Pose & pose ) {
		core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( pose ));
		task->restrict_to_repacking();
		return task;
	}

	///@brief check what happens to a packer task when the fix_his function is called
	core::pack::task::PackerTaskOP
	task_repack_with_fixhis( core::pose::Pose & pose ) {
		core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( pose ));
		utility::vector1<int> const positions;
		task->restrict_to_repacking().or_fix_his_tautomer(positions, true);
		return task;
	}

	///@brief check what happens to a packer task when the resfile fix_his command is used
	core::pack::task::PackerTaskOP
	task_repack_with_resfile( core::pose::Pose & pose ) {
		core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( pose ));
		task->restrict_to_repacking();
		core::pack::task::parse_resfile_string( pose, *task, "NATAA FIX_HIS_TAUTOMER \n start");
		return task;
	}

	// loop over 4 residues and check if the allowed residue types in the packer task
	// contain (or lack) the indicated proton
	void assert_restype_w_atom_presence( core::pack::task::PackerTaskOP task, std::string const & proton, bool presence )
	{
		for(core::Size ii = 2; ii < task->total_residue(); ++ii ) {
			bool present = false;
			for ( core::pack::task::ResidueLevelTask::ResidueTypeCOPListConstIter
					allowed_iter = task->residue_task( ii ).allowed_residue_types_begin(),
					allowed_end = task->residue_task( ii ).allowed_residue_types_end();
					allowed_iter != allowed_end; ++allowed_iter ) {
				if ( (*allowed_iter)->aa() != core::chemical::aa_his ) continue;
				if ( (*allowed_iter)->has(proton) ) {
					present = true;
				}
			}
			TS_ASSERT_EQUALS( present, presence);
		}
	}

	// --------------- Test Cases --------------- //

	void test_fix_his_tautomer() {
		//// APL TEMP
		//for ( core::Size ii = 1; ii <= 100; ++ii ) {
		//	std::cout << "test_fix_his_taut " << ii << std::endl;
		//	setUp();


		core::pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, "core/pack/task/short_his.pdb" );
		core::pose::Pose const poseconst(pose);

		//debugging
		//pose.dump_pdb("unmodified.pdb");

		//get typeset and HIS and HIS_D
		core::chemical::ResidueTypeSet & typeset(core::chemical::ChemicalManager::get_instance()->nonconst_residue_type_set(core::chemical::FA_STANDARD));
		core::chemical::ResidueType const & HIS(typeset.name_map("HIS"));
		core::chemical::ResidueType const & HIS_D(typeset.name_map("HIS_D"));
		utility::vector1< core::conformation::ResidueOP > HISs, HIS_Ds;

		//create concrete HIS(_D) residues
		for(core::Size i=2; i<pose.total_residue(); ++i) HISs.push_back(core::conformation::ResidueFactory::create_residue(HIS, pose.residue(i), pose.conformation(), true));
		for(core::Size i=2; i<pose.total_residue(); ++i) HIS_Ds.push_back(core::conformation::ResidueFactory::create_residue(HIS_D, pose.residue(i), pose.conformation(), true));

		//force all nonterminal to HIS
		core::pose::Pose HISpose(poseconst);
		for(core::Size i=2; i<pose.total_residue(); ++i) HISpose.replace_residue(i, *(HISs[i-1]), true);
		//HISpose.dump_pdb("HISpose.pdb");

		//force all nonterminal to HIS_D
		core::pose::Pose HIS_Dpose(poseconst);
		for(core::Size i=2; i<pose.total_residue(); ++i) HIS_Dpose.replace_residue(i, *(HIS_Ds[i-1]), true);
		//HIS_Dpose.dump_pdb("HIS_Dpose.pdb");

		//make copies of poses for independent manipulation
		core::pose::Pose HISpose_r(HISpose);
		core::pose::Pose HISpose_rf(HISpose);
		core::pose::Pose HISpose_rr(HISpose);

		core::pose::Pose HIS_Dpose_r(HIS_Dpose);
		core::pose::Pose HIS_Dpose_rf(HIS_Dpose);
		core::pose::Pose HIS_Dpose_rr(HIS_Dpose);

		//allow them to be repacked under appropriate conditions
		core::pack::task::PackerTaskOP task_HISpose_r    = task_repack( HISpose_r );
		core::pack::task::PackerTaskOP task_HIS_Dpose_r  = task_repack( HIS_Dpose_r );
		core::pack::task::PackerTaskOP task_HISpose_rf   = task_repack_with_fixhis( HISpose_rf );
		core::pack::task::PackerTaskOP task_HIS_Dpose_rf = task_repack_with_fixhis( HIS_Dpose_rf );
		core::pack::task::PackerTaskOP task_HISpose_rr   = task_repack_with_resfile( HISpose_rr );
		core::pack::task::PackerTaskOP task_HIS_Dpose_rr = task_repack_with_resfile( HIS_Dpose_rr );

		//check for proton presence in fixed cases
		assert_restype_w_atom_presence( task_HISpose_rf, " HD1", false);
		assert_restype_w_atom_presence( task_HISpose_rf, " HE2", true);
		assert_restype_w_atom_presence( task_HISpose_rr, " HD1", false);
		assert_restype_w_atom_presence( task_HISpose_rr, " HE2", true);
		assert_restype_w_atom_presence( task_HIS_Dpose_rf, " HE2", false);
		assert_restype_w_atom_presence( task_HIS_Dpose_rf, " HD1", true);
		assert_restype_w_atom_presence( task_HIS_Dpose_rr, " HE2", false);
		assert_restype_w_atom_presence( task_HIS_Dpose_rr, " HD1", true);

		//check for variable tautomer in nonfixed cases
		assert_restype_w_atom_presence( task_HISpose_r, " HD1", true);
		assert_restype_w_atom_presence( task_HISpose_r, " HE2", true);
		assert_restype_w_atom_presence( task_HIS_Dpose_r, " HD1", true);
		assert_restype_w_atom_presence( task_HIS_Dpose_r, " HE2", true);

	}

};//end class
