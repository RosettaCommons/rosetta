// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/task_operations/CreateSequenceMotifOperationTests.cxxtest.hh
/// @brief  Test suite for motif operation.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/task_operations/SequenceMotifTaskOperation.hh>
#include <core/select/residue_selector/ReturnResidueSubsetSelector.hh>
// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/TaskFactory.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>

static basic::Tracer TR("CreateSequenceMotifOperationTests");

using namespace core::pack::task;
using namespace protocols::task_operations;
using namespace core::select::residue_selector;

class CreateSequenceMotifOperationTests : public CxxTest::TestSuite {
	//Define Variables



public:

	core::pose::Pose pose_;
	std::string UT_inpath;
	std::string UT_outpath;

	void setUp(){
		core_init();
		core::import_pose::pose_from_file(pose_, "protocols/antibody/aho_with_antigen.pdb", core::import_pose::PDB_file);

		UT_inpath = "protocols/task_operations/";
		UT_outpath = "/Users/jadolfbr/Documents/Rosetta/main/source/test/ut_files";

	}

	void tearDown(){

	}
	void test_task_op(){

		utility::vector1< bool > subset(pose_.size(), false);

		subset[160] = true; //24L, start of L1

		ReturnResidueSubsetSelectorOP selector = ReturnResidueSubsetSelectorOP( new ::ReturnResidueSubsetSelector( subset ));

		std::string motif = "VT[^S][%PROPERTY AROMATIC][ST]-X";
		TaskFactoryOP tf = TaskFactoryOP( new TaskFactory() );
		SequenceMotifTaskOperationOP motif_op = SequenceMotifTaskOperationOP( new SequenceMotifTaskOperation( selector, motif ) );


		tf->push_back(motif_op);
		output_or_test(tf, pose_, false, "SequenceMotifOperation");


	}

	///@brief Output a UT file (pass first run) or check that the task output matches the UT file.
	void output_or_test(core::pack::task::TaskFactoryCOP tf,core::pose::Pose const & pose, bool first_run, std::string name) {


		std::string inname = UT_inpath+"/"+name+".u";
		std::string outname = UT_outpath+"/"+name+".u";
		if ( first_run ) {
			TR <<"////"<<  std::endl << inname << std::endl << std::endl;
			PackerTaskOP task = tf->create_task_and_apply_taskoperations(pose);
			task->show(TR);

			//Write to a file, which we can manually check, diff, and rename if needed.
			std::ofstream OUT;

			TR << outname << std::endl;
			OUT.open(outname.c_str());
			task->show(OUT);
			OUT << std::endl;
			OUT.close();

		} else {
			test::UTracer UT(inname);
			tf->create_task_and_apply_taskoperations(pose)->show(UT);
		}
	}








};
