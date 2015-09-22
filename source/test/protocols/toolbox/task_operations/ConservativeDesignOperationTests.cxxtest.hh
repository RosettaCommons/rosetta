// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/antibody/constraints/AntibodyConstraintTests.cxxtest.hh
/// @brief Test ConservativeDesignOperation
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

//#define private public
//#define protected public

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/util.hh>

#include <protocols/toolbox/task_operations/ConservativeDesignOperation.hh>
#include<protocols/antibody/design/util.hh>


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

#include <iostream>

#define BFE BOOST_FOREACH

using namespace protocols::antibody;
using namespace protocols::antibody::design;
using namespace core::pack::task::operation;
using namespace core::pack::task;
using namespace protocols::toolbox::task_operations;

using utility::vector1;

static THREAD_LOCAL basic::Tracer TR("protocols.antibody.ConservativeDesignOperationTests");

class ConservativeDesignOperationTests: public CxxTest::TestSuite {

public:

	core::pose::Pose pose;
	AntibodyInfoOP ab_info;
	ConservativeDesignOperationOP conservative_design;
	ConservativeDesignOperationOP conservative_design_blosum62;

	utility::vector1<core::Size> positions; //Limit to L1

	bool first_run;
	std::string UT_inpath;

	//Outpath for new UTracers. Sorry for the full path -
	// not sure a better way to do this other than copying each .u from the std output
	std::string UT_outpath;

public:

	void setUp(){
		core_init();
		core::import_pose::pose_from_pdb(pose, "protocols/antibody/aho_with_antigen.pdb"); //AHO renumbered pose
		ab_info = AntibodyInfoOP( new AntibodyInfo(pose, AHO_Scheme, North) );
		conservative_design = ConservativeDesignOperationOP( new ConservativeDesignOperation());
		conservative_design->set_data_source("chothia_1976");
		conservative_design_blosum62 = ConservativeDesignOperationOP( new ConservativeDesignOperation("BLOSUM62"));

		//Setup positions.
		for ( core::Size i = ab_info->get_CDR_start(l1, pose); i <= ab_info->get_CDR_end(l1, pose); ++i ) {
			positions.push_back(i);
		}
		conservative_design->limit_to_positions(positions);
		conservative_design_blosum62->limit_to_positions(positions);

		//Test output
		first_run = false;
		UT_inpath = "protocols/toolbox/task_operations/";
		UT_outpath = "/Users/jadolfbr/Documents/modeling/rosetta/Rosetta/main/source/test/ut_files";



	}

	void test_task_op() {



		TaskFactoryOP tf( new TaskFactory());

		//Only allow L1 to design for simplicity.
		// L1 will use conservative design
		// H3 will design into anything.

		RestrictResidueToRepackingOP restrict_l2 = protocols::antibody::design::disable_design_cdr(ab_info, l2, pose);
		RestrictResidueToRepackingOP restrict_l3 = protocols::antibody::design::disable_design_cdr(ab_info, l3, pose);
		RestrictResidueToRepackingOP restrict_h1 = protocols::antibody::design::disable_design_cdr(ab_info, h1, pose);
		RestrictResidueToRepackingOP restrict_h2 = protocols::antibody::design::disable_design_cdr(ab_info, h2, pose);
		RestrictResidueToRepackingOP restrict_h3 = protocols::antibody::design::disable_design_cdr(ab_info, h3, pose);

		RestrictResidueToRepackingOP restrict_antigen = protocols::antibody::design::disable_design_antigen(ab_info, pose);
		RestrictResidueToRepackingOP restrict_framework = protocols::antibody::design::disable_design_framework(ab_info, pose);

		tf->push_back(restrict_l2);
		tf->push_back(restrict_l3);
		tf->push_back(restrict_h1);
		tf->push_back(restrict_h2);
		tf->push_back(restrict_h3);
		tf->push_back(restrict_antigen);
		tf->push_back(restrict_framework);

		TaskFactoryOP all_but_cons = tf->clone();

		//Test Replacing, with natives.
		tf->push_back(conservative_design);
		output_or_test(tf, pose, first_run, "ConservativeDesignOperation_defaults");


		//Test Replacing without natives.
		conservative_design->include_native_aa(false);
		tf =all_but_cons->clone();
		tf->push_back(conservative_design);
		output_or_test(tf, pose, first_run, "ConservativeDesignOperation_no_native");


		//Test using a different data source as the conservative mutational data.
		tf = all_but_cons->clone();
		tf->push_back(conservative_design_blosum62);
		output_or_test(tf, pose, first_run, "ConservativeDesignOperation_blosum62");

	}

	///@brief Output a UT file (pass first run) or check that the task output matches the UT file.
	void output_or_test(core::pack::task::TaskFactoryCOP tf,core::pose::Pose const & pose, bool first_run, std::string name) {


		std::string inname = UT_inpath+"/"+name+".u";
		std::string outname = UT_outpath+"/"+name+".u";
		if ( first_run ) {
			TR <<"////"<<  std::endl << inname << std::endl << std::endl;
			PackerTaskOP task = tf->create_task_and_apply_taskoperations(pose);
			task->show(std::cout);

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

