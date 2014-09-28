// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/antibody/constraints/AntibodyConstraintTests.cxxtest.hh
/// @brief Test creation of TF and TaskOps by SeqDesignTFCreator
/// @author Jared Adolf-Bryfogle

//#define private public
//#define protected public

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/clusters/CDRClusterEnum.hh>

#include <protocols/antibody/design/CDRGraftDesignOptions.hh>
#include <protocols/antibody/design/CDRSeqDesignOptions.hh>
#include <protocols/antibody/database/CDRSetOptions.hh>
#include <protocols/antibody/database/CDRSetOptionsParser.hh>
#include <protocols/antibody/database/AntibodyDatabaseManager.hh>
#include <protocols/antibody/design/AntibodySeqDesignTFCreator.hh>
#include <protocols/antibody/design/util.hh>

#include <protocols/toolbox/task_operations/RestrictToLoopsAndNeighbors.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <basic/Tracer.hh>
#include <boost/foreach.hpp>
#include <utility/vector1.hh>

#define BFE BOOST_FOREACH

using namespace protocols::antibody;
using namespace protocols::antibody::design;
using namespace core::pack::task::operation;
using namespace core::pack::task;
using namespace protocols::toolbox::task_operations;

using utility::vector1;

static thread_local basic::Tracer TR("protocols.antibody.AntibodySeqDesign");

class AntibodySeqDesign: public CxxTest::TestSuite {
	
public:
	
	core::pose::Pose pose;
	AntibodyInfoOP ab_info;
	core::scoring::ScoreFunctionOP score;
	utility::vector1<CDRSeqDesignOptionsOP> design_options;
	
public:
	
	void setUp(){
		core_init();
		core::import_pose::pose_from_pdb(pose, "protocols/antibody/aho_with_antigen.pdb"); //AHO renumbered pose
		ab_info = AntibodyInfoOP( new AntibodyInfo(pose, AHO_Scheme, North) );
		score = core::scoring::get_score_function();
		score->score(pose); //Segfault prevention from tenA neighbor graph.  Should be a better way then this...
		
		CDRSeqDesignOptionsOP l1_options( new CDRSeqDesignOptions(l1) );
		CDRSeqDesignOptionsOP l2_options( new CDRSeqDesignOptions(l2) );
		CDRSeqDesignOptionsOP l3_options( new CDRSeqDesignOptions(l3) );
		
		CDRSeqDesignOptionsOP h1_options( new CDRSeqDesignOptions(h1) );
		CDRSeqDesignOptionsOP h2_options( new CDRSeqDesignOptions(h2) );
		CDRSeqDesignOptionsOP h3_options( new CDRSeqDesignOptions(h3) );
		
		l1_options->design(true);
		l2_options->design(true);
		l3_options->design(true);
		h1_options->design(true);
		h2_options->design(false);
		h3_options->design(true);
		
		l1_options->design_strategy(seq_design_profiles);
		l2_options->design_strategy(seq_design_profiles);
		l3_options->design_strategy(seq_design_profiles);
		h1_options->design_strategy(seq_design_conservative);
		h2_options->design_strategy(seq_design_conservative);
		h3_options->design_strategy(seq_design_basic);
		
		design_options.clear();
		design_options.resize(6, NULL);
		
		design_options[ l1 ] = l1_options;
		design_options[ l2 ] = l2_options;
		design_options[ l3 ] = l3_options;
		design_options[ h1 ] = h1_options;
		design_options[ h2 ] = h2_options;
		design_options[ h3 ] = h3_options;
		
	}

	void test_task_ops() {
		
		//Manual control of options - move this to setup.  Move each test to its own function.
	
		AntibodySeqDesignTFCreator creator =  AntibodySeqDesignTFCreator(ab_info, design_options, true);
		TaskFactoryOP tf( new TaskFactory() );
		
		TR << "--CDRs with all of them designing" << std::endl;
		RestrictToLoopsAndNeighborsOP all_cdrs = creator.generate_task_op_all_cdr_design(pose, false);
		tf->push_back(all_cdrs);
		utility::vector1<bool> disabled_cdrs(6, false);
		assert_cdr_design_is_enabled_or_disabled(tf->create_task_and_apply_taskoperations(pose), disabled_cdrs);
		
		TR << "--CDRs with only those from options designing" << std::endl;
		RestrictToLoopsAndNeighborsOP design_cdrs = creator.generate_task_op_cdr_design(pose, false);
		tf->clear();
		tf->push_back(design_cdrs);
		disabled_cdrs[ h2 ] = true;
		assert_cdr_design_is_enabled_or_disabled(tf->create_task_and_apply_taskoperations(pose), disabled_cdrs);
		
		TR<< "--CDRs designing based on vector1 bool." << std::endl;
		utility::vector1<bool> designing(6, true);
		designing[ 1 ] = false;
		RestrictToLoopsAndNeighborsOP design_vec_cdrs = creator.generate_task_op_cdr_design(pose, designing, false);
		tf->clear();
		tf->push_back(design_vec_cdrs);
		designing.flip();
		assert_cdr_design_is_enabled_or_disabled(tf->create_task_and_apply_taskoperations(pose), designing);
		
		
	
	}
	void test_normal_tf_generation() {
		AntibodySeqDesignTFCreator creator =  AntibodySeqDesignTFCreator(ab_info, design_options, true);
		TaskFactoryOP tf( new TaskFactory() );
		
		TR << "--General TF - no limts on packing or design" << std::endl;
		creator.design_antigen(false);
		creator.design_framework(false);
		utility::vector1<bool> disabled_cdrs(6, false);
		disabled_cdrs[ h2 ] = true;
		
		tf = creator.generate_tf_seq_design(pose);
		PackerTaskOP task = tf->create_task_and_apply_taskoperations(pose);
		assert_region_design_is_disabled(task, antigen_region);
		assert_region_design_is_disabled(task, framework_region);
		
		///Make sure conservative and profile-based sequences match up.
		///Extremely complicated without the diff.  
		task->show(std::cout);
		
		TR << "--TF disable antigen neighbors" << std::endl;
		creator.design_antigen(true);
		tf = creator.generate_tf_seq_design(pose);
		assert_region_design_is_disabled(tf->create_task_and_apply_taskoperations(pose), framework_region);
		
		
		TR << "--TF disable framework neighbors" << std::endl;
		creator.design_framework(true);
		creator.design_antigen(false);
		tf = creator.generate_tf_seq_design(pose);
		assert_region_design_is_disabled(tf->create_task_and_apply_taskoperations(pose), antigen_region);
		
		
	}
	void test_graft_tf_generation() {
		AntibodySeqDesignTFCreator creator =  AntibodySeqDesignTFCreator(ab_info, design_options, true);
		TaskFactoryOP tf( new TaskFactory() );
		
		TR << "--TF used for graft design" << std::endl;
		utility::vector1<bool> h1_neighbors(6, false);
		h1_neighbors[ h2 ] = true;
		h1_neighbors[ h3 ] = true;
		
		creator.design_framework(true);
		creator.design_antigen(false);
		
		tf = creator.generate_tf_seq_design_graft_design(pose, h1, h1_neighbors);
		assert_region_design_is_disabled(tf->create_task_and_apply_taskoperations(pose), antigen_region);
		
		utility::vector1<bool> disabled_cdrs(6, false);
		disabled_cdrs[ h2 ] = true; //Goes with options above.
		assert_cdr_design_disabled(tf->create_task_and_apply_taskoperations(pose), disabled_cdrs);
		
	}
	void test_utility_functions() {
		
		RestrictResidueToRepackingOP disable_antigen = disable_design_region(ab_info, pose, antigen_region);
		RestrictResidueToRepackingOP disable_framework = disable_design_region(ab_info, pose, framework_region);
		RestrictResidueToRepackingOP disable_cdrs = disable_design_region(ab_info, pose, cdr_region);
		
		assert_region_design_is_disabled(disable_antigen, antigen_region);
		assert_region_design_is_disabled(disable_framework, framework_region);
		assert_region_design_is_disabled(disable_cdrs, cdr_region);
		
		RestrictResidueToRepackingOP disable_antigen2 = disable_design_antigen(ab_info, pose);
		RestrictResidueToRepackingOP disable_framework2 = disable_design_framework(ab_info, pose);
		RestrictResidueToRepackingOP disable_cdrs2 = disable_design_cdrs(ab_info, pose);
		
		assert_region_design_is_disabled(disable_antigen2, antigen_region);
		assert_region_design_is_disabled(disable_framework2, framework_region);
		assert_region_design_is_disabled(disable_cdrs2, cdr_region);

		
	}
	
	///Utility functions
	void assert_region_design_is_disabled(RestrictResidueToRepackingOP disable, AntibodyRegionEnum region){
		TaskFactoryOP tf( new TaskFactory() );
		tf->push_back(disable);
		assert_region_design_is_disabled(tf->create_task_and_apply_taskoperations(pose), region);
	}
	void assert_region_design_is_disabled(PackerTaskOP task, AntibodyRegionEnum region) {
		
		std::string r;
		if (region == antigen_region) {r = "antigen";};
		if (region == framework_region) {r = "framework";}
		if (region == cdr_region) {r = "cdr";}
		
		TR <<"Checking region: " << r << std::endl;
		for (core::Size i = 1; i <= pose.total_residue(); ++i ){
			if (ab_info->get_region_of_residue(pose, i) == region){
				TS_ASSERT_EQUALS(task->design_residue( i ), false );
			}
		}
	}
	void assert_cdr_is_disabled(RestrictResidueToRepackingOP disable, utility::vector1<bool> cdrs_to_check){
		TaskFactoryOP tf( new  TaskFactory() );
		tf->push_back(disable);
		assert_cdr_design_is_enabled_or_disabled(tf->create_task_and_apply_taskoperations(pose), cdrs_to_check);
		
	}
	void assert_cdr_design_is_enabled_or_disabled(
		PackerTaskOP task,
		utility::vector1<bool> cdrs_to_check_disabled )
	{
		//Checks to make sure that the cdrs that are set to be disabled are disabled and that the other cdr residues are ALL enabled. 
		TR << "Checking CDR Design" << std::endl;
		assert(cdrs_to_check_disabled.size() == 6);

		for (core::Size i = 1; i <= 6; ++i ){
			CDRNameEnum cdr = static_cast<CDRNameEnum>( i );
			core::Size start = ab_info->get_CDR_start(cdr, pose);
			core::Size end = ab_info->get_CDR_end(cdr, pose);
			TR <<"CDR: " << ab_info->get_CDR_name(cdr) << std::endl;
			
			for (core::Size res = start; res <= end; ++res ){
				if (cdrs_to_check_disabled [ i ]) {
					TS_ASSERT_EQUALS(task->design_residue( res ) , false);
				}
				else {
					TS_ASSERT_EQUALS(task->design_residue( res ), true);
				}	
			}
		}
		
	}
	void assert_cdr_design_disabled(
		PackerTaskOP task,
		utility::vector1<bool> cdrs_to_check_disabled)
	{
		assert(cdrs_to_check_disabled.size() == 6);
		TR << "Checking CDR Design" << std::endl;
		for (core::Size i = 1; i <= 6; ++i ){
			CDRNameEnum cdr = static_cast<CDRNameEnum>( i );
			core::Size start = ab_info->get_CDR_start(cdr, pose);
			core::Size end = ab_info->get_CDR_end(cdr, pose);
			for (core::Size res = start; res <= end; ++res ){
				if (cdrs_to_check_disabled [ i ]) {
					TS_ASSERT_EQUALS(task->design_residue( res ) , false);
				}
			}
		}
	}
	
	
	void assert_cdr_packing_is_enabled_or_disabled(
		PackerTaskOP task,
		utility::vector1<bool> cdrs_to_check_disabled)
	{
		assert(cdrs_to_check_disabled.size() == 6);
		TR << "Checking packing " << std::endl;
		for (core::Size i = 1; i <= 6; ++i ){
			CDRNameEnum cdr = static_cast<CDRNameEnum>( i );
			core::Size start = ab_info->get_CDR_start(cdr, pose);
			core::Size end = ab_info->get_CDR_end(cdr, pose);
			for (core::Size res = start; res <= end; ++res ){
				if (cdrs_to_check_disabled [ i ]) {
					TS_ASSERT_EQUALS(task->pack_residue( res ) , false);
				}
				else {
					TS_ASSERT_EQUALS(task->pack_residue( res ), true);
				}	
			}
		}
		
	}
};
