// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       apps/pilot/ralford/mp_sfxn_integration.cc
///
/// @brief      Report Features of Residues 
/// @details    Report neighbor counts for 6, 10, and 12A radii, SASA scoring for 1.4A radii, 
///				residue layer assignemnts for 5 layer scheme, residue layer assignemnts for 2 layer
///				scheme
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (4/15/14)

// App headers
#include <devel/init.hh> 

// Options System
#include <basic/options/option.hh> 
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <protocols/membrane/MembraneUnitTestMover.hh>

// Package Headers
#include <core/pose/Pose.hh> 
#include <core/types.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh> 
#include <core/scoring/ScoreFunctionFactory.hh> 

#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh> 

// Basic headers
#include <basic/Tracer.hh> 
#include <numeric/xyzVector.hh> 

// C++ Headers
#include <cstdlib>
#include <algorithm> 

static basic::Tracer TR("mp_sfxn_integration"); 

/// @brief Main Function
int main( int argc, char* argv[] )
{
    try {
        
        using namespace core::scoring; 
        using namespace protocols::membrane;
        using namespace protocols::jd2;

        using namespace basic::options;

        // Initialize Options System, RG, and All Factory_Registrators
        devel::init(argc, argv);
        
        TR << "Membrane Energy Funciton Initial Testing" << std::endl;
        TR << "@ralford 4/21/14" << std::endl;
        TR << "==============================================================" << std::endl;
		
        TR << "Setting up Membrane Protein from Initial Resources..." << std::endl;

        // Initialize Membrane protein From Membrane Mover
        MembraneUnitTestMoverOP mp = new MembraneUnitTestMover();
        mp->register_options();
        mp->init_from_cmd();
        JobDistributor::get_instance()->go(mp);
        core::pose::PoseOP membrane_pose = mp->get_membrane_pose();

        core::Vector center(-4.494, -0.776, 0.302);
        core::Vector normal(-0.240, 0.209, -0.948);

        membrane_pose->conformation().membrane_info()->set_whole_pose_embedding(center, normal);

        // Create new scoring function from membrane score weights
        ScoreFunctionOP sfxn = new ScoreFunction();

        sfxn->set_weight(core::scoring::MPEnv, 1.0);
        sfxn->set_weight(core::scoring::Menv, 1.0);

        sfxn->set_weight(core::scoring::MPPair, 1.0);
        sfxn->set_weight(core::scoring::Mpair, 1.0);

        sfxn->set_weight(core::scoring::MPCbeta, 1.0);
        sfxn->set_weight(core::scoring::Mcbeta, 1.0);

        sfxn->set_weight( MPTermini, 1.0);
        sfxn->set_weight( Menv_termini, 1.0);

        sfxn->set_weight( MPNonHelix, 1.0);
        sfxn->set_weight( Menv_non_helix, 1.0);

        sfxn->set_weight( MPTMProj, 1.0);
        sfxn->set_weight( Menv_tm_proj, 1.0);        

        sfxn->score(*membrane_pose);

        // Show Pose energies 
        membrane_pose->energies().show(std::cout);
        
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cout << "caught exception " << e.msg() << std::endl;
    }

}
