// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/antibody/constraints/AntibodyConstraintTests.cxxtest.hh
/// @brief  tests for the Antibody Design options classes and parsers
/// @author Jared Adolf-Bryfogle


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/protocols/antibody/utilities.hh>

// Project Headers
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/clusters/CDRClusterEnum.hh>

#include <protocols/antibody/task_operations/DisableCDRsOperation.hh>
#include <protocols/antibody/task_operations/DisableAntibodyRegionOperation.hh>
#include <protocols/antibody/task_operations/RestrictToCDRsAndNeighbors.hh>
#include <protocols/antibody/task_operations/AddCDRProfilesOperation.hh>
#include <protocols/antibody/task_operations/AddCDRProfileSetsOperation.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Protocol Headers
#include <basic/Tracer.hh>
#include <boost/foreach.hpp>

#include <utility/vector1.hh>

#define BFE BOOST_FOREACH

using namespace protocols::antibody;
using namespace protocols::antibody::clusters;
using namespace protocols::antibody::task_operations;
using namespace protocols::antibody::design;
using namespace core::pack::task::operation;
using namespace core::pack::task;
using utility::vector1;

static THREAD_LOCAL basic::Tracer TR("protocols.antibody.task_operations.AntibodyTaskOps");

class AntibodyTaskOps: public CxxTest::TestSuite {

	core::pose::Pose pose_;
	core::pose::Pose pose_chothia_;

	AntibodyInfoOP ab_info_;
	AntibodyInfoOP ab_info_chothia_;
	bool first_run_;

	//Outpath for new UTracers.  Can't figure out a better way then manually setting
	// it here or going deep in the compiled code as pwd for test.py is not /source
	// Please email me if you know a better way. -JAB.
	std::string first_run_outpath_;

	std::string inpath_;

public:

	void setUp(){
		core_init();
		core::import_pose::pose_from_file(pose_, "protocols/antibody/aho_with_antigen.pdb", core::import_pose::PDB_file); //AHO renumbered pose
		core::import_pose::pose_from_file(pose_chothia_, "protocols/antibody/1bln_AB_chothia.pdb", core::import_pose::PDB_file);

		//TenA NeighborGraph setup
		core::scoring::ScoreFunctionOP score = core::scoring::get_score_function();
		score->score(pose_);
		score->score(pose_chothia_);

		inpath_ = "protocols/antibody/task_operations";
		first_run_outpath_ = "/home/jadolfbr/Documents/modeling/rosetta/Rosetta/main/source/test/ut_files";
		first_run_ = false;
		ab_info_ = AntibodyInfoOP( new AntibodyInfo(pose_, AHO_Scheme, North) );
		ab_info_chothia_ = AntibodyInfoOP( new AntibodyInfo(pose_chothia_));

	}
	void test_DisableCDRsOperation(){

		utility::vector1< bool > cdrs(8, false);
		utility::vector1< bool > all_cdrs(8, true);

		cdrs[ l1 ] = true;
		cdrs[ l2 ] = true;

		all_cdrs[ l4 ] = false;
		all_cdrs[ h4 ] = false;

		//Test Constructions

		//This won't work due to Aho scheme antibody - so we need to use one that matches the default numbering scheme.

		DisableCDRsOperationOP disable_cdrs_op_default = DisableCDRsOperationOP(new DisableCDRsOperation());

		TaskFactoryOP task = TaskFactoryOP(new TaskFactory());
		task->push_back(disable_cdrs_op_default);
		assert_cdr_packing_is_enabled_or_disabled(pose_chothia_, task->create_task_and_apply_taskoperations(pose_chothia_), ab_info_chothia_, all_cdrs);

		DisableCDRsOperationOP disable_cdrs_op_ab_info = DisableCDRsOperationOP(new DisableCDRsOperation(ab_info_));
		disable_cdrs_op_ab_info->set_cdrs(cdrs);

		DisableCDRsOperationOP disable_cdrs_op_cdrs = DisableCDRsOperationOP(new DisableCDRsOperation(ab_info_, cdrs));

		utility::vector1<DisableCDRsOperationOP> default_ops;
		default_ops.push_back(disable_cdrs_op_ab_info);
		default_ops.push_back(disable_cdrs_op_cdrs);

		for ( core::Size i = 1; i <= default_ops.size(); ++i ) {
			TaskFactoryOP task2 = TaskFactoryOP(new TaskFactory);
			task2->push_back(default_ops[ i ]);
			assert_cdr_packing_is_enabled_or_disabled(pose_, task2->create_task_and_apply_taskoperations(pose_), ab_info_, cdrs);

		}


		DisableCDRsOperationOP disable_cdrs_op_only_pack = DisableCDRsOperationOP(new DisableCDRsOperation(ab_info_, cdrs, false));
		DisableCDRsOperationOP disable_cdrs_op_only_pack2 = DisableCDRsOperationOP(new DisableCDRsOperation(ab_info_, cdrs));
		disable_cdrs_op_only_pack2->set_disable_packing_and_design(false);

		utility::vector1<DisableCDRsOperationOP> design_ops;
		default_ops.push_back(disable_cdrs_op_only_pack);
		default_ops.push_back(disable_cdrs_op_only_pack2);

		for ( core::Size i = 1; i <= default_ops.size(); ++i ) {
			TaskFactoryOP task2 = TaskFactoryOP(new TaskFactory);
			task2->push_back(default_ops[ i ]);
			assert_cdr_design_is_enabled_or_disabled(pose_, task2->create_task_and_apply_taskoperations(pose_), ab_info_, cdrs);

		}

	}
	void test_DisableAntibodyRegionOperation(){


		DisableAntibodyRegionOperationOP default_op = DisableAntibodyRegionOperationOP( new DisableAntibodyRegionOperation());
		TaskFactoryOP task = TaskFactoryOP( new TaskFactory());
		task->push_back(default_op);
		assert_region_packing_is_disabled(pose_chothia_, task->create_task_and_apply_taskoperations(pose_chothia_), ab_info_chothia_, cdr_region);

		DisableAntibodyRegionOperationOP default_ab_info = DisableAntibodyRegionOperationOP( new DisableAntibodyRegionOperation(ab_info_));
		DisableAntibodyRegionOperationOP default_region = DisableAntibodyRegionOperationOP( new DisableAntibodyRegionOperation(ab_info_, cdr_region));

		utility::vector1<DisableAntibodyRegionOperationOP> default_ops;
		default_ops.push_back(default_ab_info);
		default_ops.push_back(default_region);

		for ( core::Size i = 1; i <= default_ops.size(); ++i ) {
			TaskFactoryOP task2 = TaskFactoryOP( new TaskFactory());
			task2->push_back(default_ops[ i ]);
			assert_region_packing_is_disabled(pose_, task2->create_task_and_apply_taskoperations(pose_), ab_info_, cdr_region);
		}

		//Test disabling only design.
		DisableAntibodyRegionOperationOP default_pack_only= DisableAntibodyRegionOperationOP( new DisableAntibodyRegionOperation(ab_info_, cdr_region, false));
		DisableAntibodyRegionOperationOP default_pack_only2 = DisableAntibodyRegionOperationOP( new DisableAntibodyRegionOperation(ab_info_, cdr_region));
		default_pack_only2->set_disable_packing_and_design(false);

		utility::vector1<DisableAntibodyRegionOperationOP> default_designs;
		default_designs.push_back(default_pack_only);
		default_designs.push_back(default_pack_only2);
		for ( core::Size i = 1; i <= default_designs.size(); ++i ) {
			TaskFactoryOP task2 = TaskFactoryOP( new TaskFactory());
			task2->push_back( default_designs[ i ]);
			assert_region_design_is_disabled(pose_, task2->create_task_and_apply_taskoperations(pose_), ab_info_, cdr_region);
		}

		//Test Regions
		default_pack_only->set_region( framework_region);
		task->clear();
		task->push_back( default_pack_only );
		assert_region_design_is_disabled(pose_, task->create_task_and_apply_taskoperations(pose_), ab_info_, framework_region);

		default_pack_only->set_region( antigen_region );
		task->clear();
		task->push_back( default_pack_only);
		assert_region_design_is_disabled(pose_, task->create_task_and_apply_taskoperations(pose_), ab_info_, antigen_region);

	}
	void test_RestrictToCDRsAndNeighbors(){

		utility::vector1< bool > cdrs(8, false);
		utility::vector1< bool > all_cdrs(8, true);
		cdrs[ l1 ] = true;
		cdrs[ l2 ] = true;

		all_cdrs[ l4 ] = false;
		all_cdrs[ h4 ] = false;

		//UTracer 1
		RestrictToCDRsAndNeighborsOP default_op = RestrictToCDRsAndNeighborsOP( new RestrictToCDRsAndNeighbors());
		TaskFactoryOP task = TaskFactoryOP( new TaskFactory());
		task->push_back(default_op);
		output_or_test(task, pose_chothia_, first_run_, "RestrictCDRsOperation_UTracer1",  inpath_, first_run_outpath_);

		//UTracer2
		RestrictToCDRsAndNeighborsOP default_ab_info = RestrictToCDRsAndNeighborsOP( new RestrictToCDRsAndNeighbors(ab_info_));
		default_ab_info->set_cdrs(cdrs);
		RestrictToCDRsAndNeighborsOP default_cdrs = RestrictToCDRsAndNeighborsOP( new RestrictToCDRsAndNeighbors(ab_info_, cdrs));

		task->clear();
		task->push_back(default_ab_info);
		output_or_test(task, pose_, first_run_, "RestrictCDRsOperation_UTracer2",  inpath_, first_run_outpath_);

		//UTracer 3
		RestrictToCDRsAndNeighborsOP default_cdr_design = RestrictToCDRsAndNeighborsOP( new RestrictToCDRsAndNeighbors( ab_info_, all_cdrs, true));

		task->clear();
		task->push_back(default_cdr_design);
		output_or_test(task, pose_, first_run_, "RestrictCDRsOperation_UTracer3",  inpath_, first_run_outpath_);

		//Utracer 4
		default_cdrs->set_neighbor_distance(13.0);
		default_cdrs->set_allow_design_neighbor_antigen(true);
		default_cdrs->set_allow_design_neighbor_framework(false);
		default_cdrs->set_stem_size(2);
		default_cdrs->set_allow_design_cdr(true);

		task->clear();
		task->push_back(default_cdrs);
		output_or_test(task, pose_, first_run_, "RestrictCDRsOperation_UTracer4",  inpath_, first_run_outpath_);

		//No Neighbors
		RestrictToCDRsAndNeighborsOP no_neighbors = RestrictToCDRsAndNeighborsOP( new RestrictToCDRsAndNeighbors(ab_info_, cdrs, false));
		no_neighbors->set_neighbor_distance(0);

		task->clear();
		task->push_back(no_neighbors);
		cdrs.flip();
		assert_cdr_packing_is_enabled_or_disabled(pose_, task->create_task_and_apply_taskoperations(pose_), ab_info_, cdrs);


	}
	void test_AddCDRProfilesOperation() {
		utility::vector1< bool > cdrs(8, false);
		utility::vector1< bool > all_cdrs(8, true);
		cdrs[ l1 ] = true;
		cdrs[ l2 ] = true;

		all_cdrs[ l4 ] = false;
		all_cdrs[ h4 ] = false;

		AddCDRProfilesOperationOP default_op = AddCDRProfilesOperationOP(new AddCDRProfilesOperation());
		default_op->set_force_north_paper_db(true);

		TaskFactoryOP task = TaskFactoryOP( new TaskFactory());
		task->push_back(default_op);
		output_or_test(task, pose_chothia_, first_run_, "AddCDRProfilesOperation_UTracer1",  inpath_, first_run_outpath_);

		AddCDRProfilesOperationOP default_ab_info = AddCDRProfilesOperationOP(new AddCDRProfilesOperation(ab_info_));
		default_ab_info->set_force_north_paper_db(true);

		default_ab_info->set_cdrs(cdrs);

		AddCDRProfilesOperationOP default_cdrs = AddCDRProfilesOperationOP( new AddCDRProfilesOperation(ab_info_, cdrs));
		default_cdrs->set_force_north_paper_db(true);

		task->clear();
		task->push_back(default_ab_info);
		output_or_test(task, pose_, first_run_, "AddCDRProfilesOperation_UTracer2",  inpath_, first_run_outpath_);

		//// Fallback as none
		AddCDRProfilesOperationOP no_fallback = AddCDRProfilesOperationOP( new AddCDRProfilesOperation(ab_info_));
		no_fallback->set_force_north_paper_db(true);

		no_fallback->set_fallback_strategy(design::seq_design_none);

		task->clear();
		task->push_back(no_fallback);
		output_or_test(task, pose_, first_run_, "AddCDRProfilesOperation_UTracer3",  inpath_, first_run_outpath_);

		/// Force Fallback
		AddCDRProfilesOperationOP forced_fallback = AddCDRProfilesOperationOP( new AddCDRProfilesOperation(ab_info_));
		forced_fallback->set_force_north_paper_db(true);
		forced_fallback->set_primary_strategy(seq_design_conservative);

		task->clear();
		task->push_back(forced_fallback);
		output_or_test(task, pose_, first_run_, "AddCDRProfilesOperation_UTracer4",  inpath_, first_run_outpath_);

		/// Test Pre-load of data
		AddCDRProfilesOperationOP pre_loaded_data = AddCDRProfilesOperationOP( new AddCDRProfilesOperation(ab_info_, cdrs));
		pre_loaded_data->set_force_north_paper_db(true);
		pre_loaded_data->pre_load_data(pose_);
		task->clear();
		task->push_back(pre_loaded_data);
		output_or_test(task, pose_, first_run_, "AddCDRProfilesOperation_UTracer5",  inpath_, first_run_outpath_);

		// Test Pre-load with profile sets.
		AddCDRProfilesOperationOP pre_loaded_sets = AddCDRProfilesOperationOP( new AddCDRProfilesOperation(ab_info_));
		pre_loaded_sets->set_force_north_paper_db(true);
		pre_loaded_sets->set_primary_strategy(seq_design_profile_sets);
		pre_loaded_sets->pre_load_data(pose_);
		task->clear();
		task->push_back(pre_loaded_sets);
		output_or_test(task, pose_, first_run_, "AddCDRProfilesOperation_UTracer6",  inpath_, first_run_outpath_);

		// Test Pre-load with profile sets and profiles.
		AddCDRProfilesOperationOP pre_loaded_combined = AddCDRProfilesOperationOP( new AddCDRProfilesOperation(ab_info_));
		pre_loaded_combined->set_force_north_paper_db(true);
		pre_loaded_combined->set_primary_strategy(seq_design_profile_sets_combined);
		pre_loaded_combined->pre_load_data(pose_);
		task->clear();
		task->push_back(pre_loaded_combined);
		output_or_test(task, pose_, first_run_, "AddCDRProfilesOperation_UTracer7",  inpath_, first_run_outpath_);
	}
	void test_AddCDRProfileSetsOperation() {
		utility::vector1< bool > cdrs(8, false);
		utility::vector1< bool > all_cdrs(8, true);
		cdrs[ l2 ] = true;

		all_cdrs[ l4 ] = false;
		all_cdrs[ h4 ] = false;
		
		TaskFactoryOP task = TaskFactoryOP( new TaskFactory());

		AddCDRProfileSetsOperationOP default_ab_info = AddCDRProfileSetsOperationOP( new AddCDRProfileSetsOperation(ab_info_));
		default_ab_info->set_cdrs(cdrs);
		default_ab_info->set_force_north_paper_db( true );

		task->clear();
		task->push_back(default_ab_info);
		output_or_test(task, pose_, first_run_, "AddCDRProfileSetsOperation_UTracer2",  inpath_, first_run_outpath_);

		AddCDRProfileSetsOperationOP default_op = AddCDRProfileSetsOperationOP( new AddCDRProfileSetsOperation());
		default_op->set_picking_rounds(5); //Should then sample multiple CDRs
		default_op->set_force_north_paper_db( true );

		task->clear();
		task->push_back(default_op);
		output_or_test(task, pose_chothia_, first_run_, "AddCDRProfileSetsOperation_UTracer1",  inpath_, first_run_outpath_);

		AddCDRProfileSetsOperationOP pick_op = AddCDRProfileSetsOperationOP( new AddCDRProfileSetsOperation(ab_info_, cdrs));
		pick_op->set_force_north_paper_db( true );
		pick_op->set_cutoff(10);
		pick_op->set_picking_rounds( 5 );

		task->clear();
		task->push_back(pick_op);
		output_or_test(task, pose_, first_run_, "AddCDRProfileSetsOperation_UTracer3",  inpath_, first_run_outpath_);

		AddCDRProfileSetsOperationOP length_op = AddCDRProfileSetsOperationOP( new AddCDRProfileSetsOperation(ab_info_, cdrs, true));
		length_op->set_force_north_paper_db( true );
		length_op->set_include_native_type( false );
		length_op->set_limit_only_to_length( true );
		length_op->set_picking_rounds( 5 );

		task->clear();
		task->push_back(length_op);
		output_or_test(task, pose_, first_run_, "AddCDRProfileSetsOperation_UTracer4",  inpath_, first_run_outpath_);

	}

};
