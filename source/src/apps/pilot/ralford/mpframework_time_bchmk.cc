// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       apps/pilot/ralford/frmwk_time_benchmark.cc
///
/// @brief      Quick time comparison between old and new framework
/// @details    Construct a standard pose and score 10x vs construct a membrane framework pose and score 10x
///             checking timing comparison + overhead.
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (2/14/14)

// App Headers
#include <devel/init.hh>

// Project Headers
#include <protocols/membrane/MembraneUnitTestMover.hh>
#include <core/membrane/MembraneProteinFactory.hh>

#include <core/membrane/util/Exceptions.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>
#include <ctime>

using basic::Error;
using basic::Warning;

using namespace core::membrane;

static basic::Tracer TR( "apps.pilot.ralford.framework_time_benchmark" );

/// @brief Main Method
int main( int argc, char* argv[] )
{

    try {

        using namespace protocols::membrane;
        using namespace protocols::jd2;
        using namespace core::scoring;

        // Initialize Options System, RG, and All Factory_Registrators
        devel::init(argc, argv);

        // Print App Info
        TR << "Membrane Framework Time Benchmark" << std::endl;
        TR << "Version 1.0" << std::endl;
        TR << "Last Modified: 2/15/14" << std::endl;

        TR << "Setting up membrane scoring function..." << std::endl;
        // Setup Membrane Scoring Function
        ScoreFunctionOP membrane_sfxn = ScoreFunctionFactory::create_score_function( "membrane_highres_Menv_smooth" );

        TR << "Setting up Membrane Protein from Initial Resources..." << std::endl;
        // Initialize Membrane Mover
        MembraneUnitTestMoverOP mp = new MembraneUnitTestMover();
        mp->register_options();
        mp->init_from_cmd();

        // Snd the mover through the job distributor
        JobDistributor::get_instance()->go(mp);

        // Time Membrane Overhead
        std::clock_t mp_start;
        double mp_duration;
        mp_start = std::clock();
        core::pose::PoseOP membrane_pose = mp->get_membrane_pose();
        mp_duration = ( std::clock() - mp_start ) / (double) CLOCKS_PER_SEC;

        // Time Non membrane pose overhead
        std::clock_t non_mp_start;
        double non_mp_duration;
        non_mp_start = std::clock();
        core::pose::PoseOP non_membrane_pose = mp->get_non_membrane_pose();
        non_mp_duration = ( std::clock() - non_mp_start ) / (double) CLOCKS_PER_SEC;

        // Time Non Membrane Pose Scoring
        std::clock_t non_mp_score_start;
        double non_mp_score_duration;
        non_mp_score_start = std::clock();
        membrane_sfxn->score(*non_membrane_pose);
        non_mp_score_duration = ( std::clock() - non_mp_score_start ) / (double) CLOCKS_PER_SEC;

        // Time Membrane Pose Scoring
        std::clock_t mp_score_start;
        double mp_score_duration;
        mp_score_start = std::clock();
        membrane_sfxn->score(*membrane_pose);
        mp_score_duration = ( std::clock() - mp_score_start ) / (double) CLOCKS_PER_SEC;



        // Show Results
        TR << "Membrane Protein Framework: Benchmark Timing Test" << std::endl;
        TR << "=======================================================================" << std::endl;
        TR << "Membrane Pose Initialization Time: " << mp_duration << std::endl;
        TR << "Membrane Pose Scoring Time (nstruct = 1):  " << mp_score_duration << std::endl;
        TR << "Non Membrane Pose Initialization Time: " << non_mp_duration << std::endl;
        TR << "Non Membrane Pose Scoring Time (nstruct =1 ): " << non_mp_score_duration << std::endl;
        TR << "Timing benchmark run complete!" << std::endl;


    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cout << "caught exception " << e.msg() << std::endl;
				return -1;
    }
}
