// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/antibody/constraints/ResidueProbDesignOperationTests.hh
/// @brief Test ResidueProbDesignOperation
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

#include <protocols/task_operations/ResidueProbDesignOperation.hh>
#include <protocols/antibody/database/AntibodyDatabaseManager.hh>
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

#define BFE BOOST_FOREACH

using namespace protocols::antibody;
using namespace protocols::antibody::design;
using namespace core::pack::task::operation;
using namespace core::pack::task;
using namespace protocols::task_operations;
using namespace core::chemical;

using utility::vector1;

static basic::Tracer TR("ResidueProbDesignOperationTests");

class ResidueProbDesignOperationTests: public CxxTest::TestSuite {
	typedef std::map< core::chemical::AA, core::Real > AAProbabilities; //Map of an amino acid and it's probability.
	typedef std::map< core::Size, AAProbabilities >ProbabilitySet; //Amino acid probabilities for a particular residue number.

public:

	core::pose::Pose pose;
	AntibodyInfoOP ab_info;
	ResidueProbDesignOperationOP prob_design;
	AntibodyDatabaseManagerOP db_manager;
	bool first_run;

public:

	void setUp(){
		core_init();
		core::import_pose::pose_from_file(pose, "protocols/antibody/aho_with_antigen.pdb", core::import_pose::PDB_file); //AHO renumbered pose
		ab_info = AntibodyInfoOP( new AntibodyInfo(pose, AHO_Scheme, North) );
		prob_design = ResidueProbDesignOperationOP( new ResidueProbDesignOperation());
		db_manager = AntibodyDatabaseManagerOP(new AntibodyDatabaseManager(ab_info, true /* force paper ab_db */));
		first_run = false; //Change to output to STD out rather then load U.
	}

	void test_task_op() {

		utility::vector1<bool> cdrs_to_use(6, false);

		cdrs_to_use[ l1 ] = true;

		ProbabilitySet prob_set;
		db_manager->load_cdr_design_data_for_cdrs(cdrs_to_use, pose, prob_set, 1);

		prob_design->set_aa_probability_set(prob_set);

		TaskFactoryOP tf( new TaskFactory());


		RestrictResidueToRepackingOP restrict_l2 = protocols::antibody::design::disable_design_cdr(ab_info, l2, pose);
		RestrictResidueToRepackingOP restrict_l3 = protocols::antibody::design::disable_design_cdr(ab_info, l3, pose);
		RestrictResidueToRepackingOP restrict_h1 = protocols::antibody::design::disable_design_cdr(ab_info, h1, pose);
		RestrictResidueToRepackingOP restrict_h2 = protocols::antibody::design::disable_design_cdr(ab_info, h2, pose);

		//Were going to allow H3 to design into anything.  Not that this is what you would actually do.

		RestrictResidueToRepackingOP restrict_antigen = protocols::antibody::design::disable_design_antigen(ab_info, pose);
		RestrictResidueToRepackingOP restrict_framework = protocols::antibody::design::disable_design_framework(ab_info, pose);

		tf->push_back(restrict_l2);
		tf->push_back(restrict_l3);
		tf->push_back(restrict_h1);
		tf->push_back(restrict_h2);
		tf->push_back(restrict_antigen);
		tf->push_back(restrict_framework);

		TaskFactoryOP all_but_prob = tf->clone();

		tf->push_back(prob_design);

		//Test Replacing, with natives.
		output_or_test(tf, pose, first_run, "protocols/task_operations/ResidueProbDesignOperation_defaults.u");

		//Test Replacing without natives.
		prob_design->set_include_native_restype(false);
		tf =all_but_prob->clone();
		tf->push_back(prob_design);
		output_or_test(tf, pose, first_run, "protocols/task_operations/ResidueProbDesignOperation_no_native.u");

		//Test Sampling Zero Probabilities
		prob_design->set_defaults();
		prob_design->set_sample_zero_probs_at(.1);

		tf =all_but_prob->clone();
		tf->push_back(prob_design);

		output_or_test(tf, pose, first_run, "protocols/task_operations/ResidueProbDesignOperation_zero_at_.05.u");

		//Test Picking rounds
		prob_design->set_defaults();
		prob_design->set_picking_rounds(5);

		tf =all_but_prob->clone();
		tf->push_back(prob_design);

		output_or_test(tf, pose, first_run, "protocols/task_operations/ResidueProbDesignOperation_picks_5.u");

		//Test Setting for one amino acid - make them all even just for the test
		AAProbabilities probs_24;
		for ( core::Size i = 1; i <= 20; ++i ) {
			core::chemical::AA a = static_cast<core::chemical::AA>(i);
			probs_24[ a ] = 1.0;
		}

		prob_design->set_defaults();
		prob_design->set_aa_probabilities( ab_info->get_CDR_start(l1, pose), probs_24);

		tf =all_but_prob->clone();
		tf->push_back(prob_design);
		output_or_test(tf, pose, first_run, "protocols/task_operations/ResidueProbDesignOperation_set_position.u");

		//Test Setting for all amino acids - make tyrosine and serine dominate, just because.
		probs_24[ aa_tyr] = 15.0;
		probs_24[ aa_ser] =15.0;

		prob_design->set_defaults();
		prob_design->set_include_native_restype(false);
		prob_design->set_overall_aa_probabilities(probs_24);

		tf =all_but_prob->clone();
		tf->push_back(prob_design);
		output_or_test(tf, pose, first_run, "protocols/task_operations/ResidueProbDesignOperation_set_all_aa_probs.u");
	}
	void output_or_test(TaskFactoryCOP tf, core::pose::Pose const & pose, bool first_run, std::string name) {
		if ( first_run ) {
			TR <<std::endl << "////"<<  std::endl << name << std::endl << std::endl;
			tf->create_task_and_apply_taskoperations(pose)->show(TR);
		} else {
			test::UTracer UT(name);
			tf->create_task_and_apply_taskoperations(pose)->show(UT);
		}
	}

};

