// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/moves/MonteCarloTests.cxxtest.hh
/// @brief  Tests for MonteCarlo and associated machinery.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/monte_carlo/MonteCarloInterface.hh>
#include <protocols/analysis/InterfaceAnalyzerMover.hh>
#include <protocols/antibody/design/GeneralAntibodyModeler.hh>
#include <protocols/antibody/AntibodyInfo.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("MonteCarloTests");

using namespace protocols::moves;
using namespace protocols::monte_carlo;
using namespace core::scoring;
using namespace protocols::antibody::design;
using namespace protocols::antibody;
using namespace protocols::analysis;

class MonteCarloTests : public CxxTest::TestSuite {
	//Define Variables

public:

	core::pose::Pose ab_pose_;
	core::pose::Pose ab_pose2_;

	ScoreFunctionOP scorefxn_;

	void setUp(){
		core_init();

		core::import_pose::pose_from_file(ab_pose_, "protocols/antibody/1bln_AB_aho.pdb", core::import_pose::PDB_file);

		core::import_pose::pose_from_file(ab_pose2_, "protocols/antibody/4LGS_AHO.pdb", core::import_pose::PDB_file);

		scorefxn_ = get_score_function();

	}

	void tearDown(){

	}



	void test_vanilla_mc(){

		//JAB - this is just testing basic functionality as none of this existed before and we are now
		// in 2018.
		MonteCarlo mc = MonteCarlo(ab_pose_, *scorefxn_, 1.0);

		TS_ASSERT( mc.last_score() != 0);
		TS_ASSERT( mc.lowest_score() != 0);
		TS_ASSERT( mc.last_accepted_score() != 0);

		bool accepted = mc.boltzmann( -10000, ab_pose2_ );
		TS_ASSERT_EQUALS( accepted, true);

		//Check coordinate equivalency to check that the last accepted pose is actually the one we set through boltzmann.

		TS_ASSERT( core::pose::compare_atom_coordinates( mc.last_accepted_pose(), ab_pose2_));
		TS_ASSERT( core::pose::compare_atom_coordinates( mc.lowest_score_pose(), ab_pose2_));

		mc.set_lowest_score_pose( ab_pose_, -11000 );

		TS_ASSERT( core::pose::compare_atom_coordinates( mc.lowest_score_pose(), ab_pose_));

		//Reset ab_pose2_ as the lowest and last accepted score pose.  Make sure it is accepted.
		core::pose::Pose bogus_pose_;
		accepted = mc.boltzmann( -20000, ab_pose2_ );
		TS_ASSERT_EQUALS( accepted, true);

		//Makes sure a bogus high score is not accepted.
		accepted = mc.boltzmann( 10000, bogus_pose_ );
		TS_ASSERT_EQUALS(accepted, false);

		//Assert that the boltzmann pose is not the lowest score pose.

		TS_ASSERT(!core::pose::compare_atom_coordinates( mc.lowest_score_pose(), ab_pose_));

		//Makes sure that the last accepted pose is still the original ab_pose.
		TS_ASSERT( core::pose::compare_atom_coordinates( mc.lowest_score_pose(), ab_pose2_));
	}

	void test_interface_mc(){

		AntibodyInfoOP ab_info = AntibodyInfoOP(new AntibodyInfo(ab_pose_, AHO_Scheme, North));
		GeneralAntibodyModeler modeler = GeneralAntibodyModeler(ab_info);
		//modeler.ab_dock_chains("L_H");
		InterfaceAnalyzerMover analyzer = InterfaceAnalyzerMover("L_H");
		analyzer.set_compute_interface_delta_hbond_unsat( false );
		analyzer.set_compute_packstat( false );
		analyzer.set_compute_interface_sc( false );
		analyzer.set_calc_hbond_sasaE( false );
		analyzer.set_use_centroid_dG( false );
		analyzer.set_pack_separated( false ); //Can't as the scores will be too different.
		analyzer.set_calc_dSASA( false );
		analyzer.set_scorefunction( scorefxn_ );
		analyzer.apply(ab_pose_);

		MonteCarloInterface mc = MonteCarloInterface(ab_pose_, *scorefxn_, 1.0, "L_H");
		mc.set_dG_weight( 1.0 );
		mc.set_repack_separated( false );
		mc.reset( ab_pose_ );

		//That called reset and set the last accepted score.  So now we make sure it's the same score
		// as interface analyzer mover.

		core::Real analyzer_score = analyzer.get_interface_dG();

		TS_ASSERT_DELTA( analyzer_score, mc.last_accepted_score(), 1.0);
		TS_ASSERT_DELTA( analyzer_score, mc.lowest_score(), 1.0);

		core::Real total_weight = 0.5;

		mc.set_total_weight( total_weight );
		mc.reset( ab_pose_);

		core::Real total_s = scorefxn_->score( ab_pose_ );

		core::Real weighted_score = (1.0 * analyzer_score) + (total_weight * total_s );
		TS_ASSERT_DELTA( weighted_score, mc.last_accepted_score(), 1.0);
		TS_ASSERT_DELTA( weighted_score, mc.lowest_score(), 1.0);

	}






};
