// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constaints/DihedralConstraint.cxxtest.hh
/// @brief  test suite for constraints between protein and ligand
/// @author Florian Richter

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
//#include <core/conformation/Residue.hh>


//#include <core/kinematics/MoveMap.hh>

//#include <core/optimization/AtomTreeMinimizer.hh>
//#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/chemical/ChemicalManager.hh> //need for additional residue
#include <core/chemical/ResidueTypeSet.hh>
#include <basic/options/option.hh> //needed to set option
#include <basic/MetricValue.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NonlocalContactsCalculator.hh>


//packing stuff
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/moves/Mover.hh>

#include <core/types.hh>


#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/run.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


using basic::Error;
using basic::Warning;

static basic::Tracer TR("HBondCalculatorsTest.cxxtest");

using namespace core;


class HBondCalculatorsTest : public CxxTest::TestSuite
{

public:
	HBondCalculatorsTest() {};


	// Shared initialization goes here.
	void setUp() {

		core_init_with_additional_options("-run:preserve_header -extra_res_fa protocols/enzdes/D2N.params");
		/*
		// Residue definitions can't be supplied on the command line b/c
		// the ResidueTypeSet is already initialized.
		using namespace core::chemical;
		utility::vector1< std::string > params_files;
		ResidueTypeSetCOP const_residue_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		ResidueTypeSet & residue_set = const_cast< ResidueTypeSet & >(*const_residue_set);
		if ( !residue_set.has_name("D2N") ) params_files.push_back("protocols/enzdes/D2N.params");
		residue_set.read_files_for_custom_residue_types(params_files);
		basic::options::option[basic::options::OptionKeys::run::preserve_header ].value(true);

		//enz_io = new protocols::enzdes::EnzConstraintIO(& residue_set);
		*/

	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_hb_calculators()
	{
		using namespace core::scoring::constraints;
		//typedef core::id::AtomID AtomID;

		core::pose::Pose test_pose;
		core::import_pose::pose_from_file( test_pose, "protocols/enzdes/ligtest_it.pdb", core::import_pose::PDB_file);
		scoring::ScoreFunctionOP scorefxn;
		scorefxn = core::scoring::get_score_function();
		(*scorefxn)(test_pose);

		time_t start1, end1, start2, end2;

		core::pose::metrics::PoseMetricCalculatorOP hbtest_calc( new protocols::toolbox::pose_metric_calculators::NumberHBondsCalculator() );

		core::pose::metrics::CalculatorFactory::Instance().register_calculator( "hbcalcname", hbtest_calc );
		core::pose::metrics::PoseMetricCalculatorOP unsattest_calc( new protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator("sasacalcname","hbcalcname") );

		core::pose::metrics::CalculatorFactory::Instance().register_calculator( "unsatcalcname", unsattest_calc );

		core::pose::metrics::PoseMetricCalculatorOP packstattest_calc( new protocols::toolbox::pose_metric_calculators::PackstatCalculator() );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( "packstatcalcname", packstattest_calc );

		core::pose::metrics::PoseMetricCalculatorOP noligpackstattest_calc( new protocols::toolbox::pose_metric_calculators::PackstatCalculator(true) );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( "noligpackstatcalcname", noligpackstattest_calc );

		core::pose::metrics::PoseMetricCalculatorOP nlcontactstest_calc( new protocols::toolbox::pose_metric_calculators::NonlocalContactsCalculator() );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( "nlcontactscalcname", nlcontactstest_calc );


		basic::MetricValue< utility::vector1< core::Size > >mval;

		time(&start1);
		test_pose.metric("hbcalcname", "residue_Hbonds", mval );
		time(&end1);

		utility::vector1< core::Size > res_hb = mval.value();

		//absolute test for number hbonds
		//TS_ASSERT_DELTA(res_hb[1], 0, 1e-3 );


		//now let's repack a few aliphatic residues
		core::pack::task::PackerTaskOP test_task;
		test_task = core::pack::task::TaskFactory::create_packer_task( test_pose );

		for ( core::Size i = 1; i <= test_pose.size(); ++i ) {
			if ( (i==22) || (i==37) || (i==59) ) test_task->nonconst_residue_task(i).restrict_to_repacking();
			else test_task->nonconst_residue_task(i).prevent_repacking();
		}

		protocols::minimization_packing::PackRotamersMoverOP testpack( new protocols::minimization_packing::PackRotamersMover(scorefxn, test_task) );
		testpack->apply(test_pose);
		(*scorefxn)(test_pose);

		time(&start2);
		test_pose.metric("hbcalcname", "residue_Hbonds", mval );
		time(&end2);

		Real dif1 = difftime(end1,start1);
		Real dif2 = difftime(end2,start2);

		TR << "first calc took "<< dif1 <<" seconds, second calc took " << dif2 << " seconds" << std::endl;

		utility::vector1< core::Size > res2_hb = mval.value();

		//after this repack step, all hbonds should still be intact
		for ( core::Size i = 1; i <= test_pose.size(); ++i ) {
			//   std::cerr << "res " << i <<", hbbef " << res_hb[i] <<", hbaft " << res2_hb[i] << std::endl;
			TS_ASSERT_DELTA( res_hb[i],res2_hb[i],1e-3);
		}


		basic::MetricValue< utility::vector1< core::Size > >res_bur_mval;
		test_pose.metric("unsatcalcname", "residue_bur_unsat_polars", res_bur_mval );

		TR << "residue 6 has " << res_bur_mval.value()[6] << " buried unsatisfied polars." << std::endl;

		//for( core::Size i = 1; i <= test_pose.size(); ++i){
		// TR << "residue " << i << " has " <<  res_bur_mval.value()[i] << " bur unsat polars." << std::endl;
		//}

		TS_ASSERT_DELTA( res_bur_mval.value()[6],2,1e-3);

		basic::MetricValue< core::Size > tot_nlcontacts_mval;
		test_pose.metric("nlcontactscalcname", "total_nlcontacts", tot_nlcontacts_mval );

		TR << "testpose has " << tot_nlcontacts_mval.value() << " non-local contacts." << std::endl;

		basic::MetricValue< core::Real > tot_packstat_mval;
		basic::MetricValue< utility::vector1< core::Real > > residue_packstat_mval;

		basic::MetricValue< core::Real > nolig_tot_packstat_mval;
		basic::MetricValue< utility::vector1< core::Real > > nolig_residue_packstat_mval;

		test_pose.metric("packstatcalcname", "total_packstat", tot_packstat_mval );
		TR << "done first packstat calc " << std::endl;
		test_pose.metric("packstatcalcname", "residue_packstat", residue_packstat_mval );
		TR << "done second packstat calc (should have been a fast lookup this time) " << std::endl;


		test_pose.metric("noligpackstatcalcname", "total_packstat", nolig_tot_packstat_mval );
		TR << "done first nolig packstat calc " << std::endl;
		test_pose.metric("noligpackstatcalcname", "residue_packstat", nolig_residue_packstat_mval );
		TR << "done second nolig packstat calc (should have been a fast lookup this time) " << std::endl;


		utility::vector1< core::Real > residue_packstat = residue_packstat_mval.value();

		core::Real respackstat_sum(0.0);
		for ( core::Size i = 1; i <= residue_packstat.size(); ++i ) {
			respackstat_sum = respackstat_sum + residue_packstat[i];
		}

		core::Real av_respackstat = respackstat_sum / residue_packstat.size();

		TR << "total packstat is " << tot_packstat_mval.value() << ", average respackstat is " << av_respackstat <<  std::endl;


		utility::vector1< core::Real > nolig_residue_packstat = nolig_residue_packstat_mval.value();

		core::Real noligrespackstat_sum(0.0);
		for ( core::Size i = 1; i <= residue_packstat.size(); ++i ) {
			noligrespackstat_sum = noligrespackstat_sum + nolig_residue_packstat[i];
		}

		core::Real nolig_av_respackstat = noligrespackstat_sum / nolig_residue_packstat.size();

		TR << "nolig total packstat is " << nolig_tot_packstat_mval.value() << ", average nolig respackstat is " << nolig_av_respackstat <<  std::endl;

	}


};

