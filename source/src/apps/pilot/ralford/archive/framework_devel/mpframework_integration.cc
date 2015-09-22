// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       apps/pilot/ralford/mpframework_integration.cc
///
/// @brief      Integration test for the Membrane Protein Framework
/// @details    The purpose of this integration test is to act as a top level test for the membrane
///             conformational framework in Rosetta. It checks conformation, fold tree construction,
///             virtual residue placement, and that no scores change by adding the conformational framework
///             otehr than due to randomization.
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (1/17/14)

/// Note - this test is not meant to be a specific unit test. It is meant to check the
/// general form and setup of membrane proteins in Rosetta. For unit tests - look in test/core/membrane

// App Headers
#include <devel/init.hh>

#include <protocols/membrane/MembraneUnitTestMover.hh>

// Project Headers
#include <core/membrane/properties/SpanningTopology.hh>

#include <core/membrane/util/definitions.hh>
#include <core/membrane/util/Exceptions.hh>

#include <core/membrane/MembraneConformation.hh>
#include <core/conformation/Conformation.hh>

#include <core/kinematics/FoldTree.hh>

#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>
#include <basic/options/keys/membrane.OptionKeys.gen.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

using basic::Error;
using basic::Warning;

using namespace core::membrane;

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.ralford.mpframework_integration" );

/// Testing Utility Functions /////////////
core::Size get_ddG( core::pose::Pose & mp_pose, core::pose::Pose & non_mp_pose, core::scoring::ScoreType type ) {

    // Grab energies
    core::Real score1 = mp_pose.energies().total_energies()[ type ];
    core::Real score2 = non_mp_pose.energies().total_energies()[ type ];

    return score1-score2;
}

///// Individual Unit Tests ////////////

/// @brief      Test Contents of the Membrane Chain
/// @details    Check that the pose contains a membrane residue
///             and an embedding data residue for every original polymer chain
///             in the pose
void test_membrane_chain( core::pose::Pose & pose ) {

    TR << "Test: Membrane Chains...";

    // Calculate the number of chains in the original pose
    core::Size nchains_orig = pose.conformation().num_chains()-1;
    core::Size nchains_final = pose.conformation().num_chains();

    // Check the contents of the last chain
    core::Size startres = pose.conformation().chain_begin(nchains_final);
    core::Size endres = pose.conformation().chain_end(nchains_final);

    // counter for nchains check
    core::Size counter = endres-startres;

    for ( core::Size i = startres; i <= endres; ++i ) {

        // Check residue
        if ( i == startres ) {
            if ( pose.residue(i).type().name().compare("MEM") != 0 ) {
                TR << "Error: First residue in membrane chain must be a membrane residue type" << std::endl;
                TR << "Error: Looking for residue type MEM, found type " << pose.residue(i).type().name() << std::endl;
                utility_exit_with_message("Integration test failed!");
            }

        } else {

            if ( pose.residue(i).type().name().compare("EMB") != 0 ) {
                TR << "Error: Non first chain residues in the membrane chain must be embedding residues" << std::endl;
                TR << "Error: Looking for residue type EMB, found type " << pose.residue(i).type().name() << std::endl;
                utility_exit_with_message("Integration test failed!");
            }
        }

        // Check total number of embedding residues == nchains orig
        if ( nchains_orig != counter ) {
            TR << "Error: Number of embedding residues is not equal to number of protein chains in the pose!" << std::endl;
            TR << "Error: Looking for " << nchains_orig << " and found " << counter << std::endl;
        }
    }

    TR << "Test Passed!" << std::endl;
}

/// @brief      Test the Contents of the Membrane Conformation
/// @details    Check that the position of the embedding residues is consistent
///             with teh numebring rules as well as the residue numbering for the membrane
///             residue. Check the contents of the embedding map is consistent with
///             membrane data. Check number of chains of topology and total res matches up
///             with the actual total residue.
void test_membrane_conformation( core::pose::Pose & pose ) {

    using namespace core::membrane;

    TR << "Test: Membrane Conformation...";

    // Grab Membrane Conformation from the pose
    MembraneConformation & mp_conf = dynamic_cast< MembraneConformation & >( pose.conformation() );

    // Pose Chain Info
    core::Size nchains_orig = pose.conformation().num_chains()-1;
    core::Size nchains_final = pose.conformation().num_chains();

    // Checking consistency in the membrane root
    if ( pose.residue( mp_conf.membrane() ).type().name().compare("MEM") ) {
        TR << "Error: Membrane root does not specify membrane residue in the pose" << std::endl;
        TR << "Error: Specified membrane root residue number " << mp_conf.membrane() << "has type name of " << pose.residue( mp_conf.membrane() ).type().name() << std::endl;
        utility_exit_with_message("Integration test failed!");
    }

    if ( (core::Size) mp_conf.membrane() != pose.conformation().chain_begin( nchains_final ) ) {
        TR << "Error: Membrane residue must be the first residue in the last chain of the membrane pose" << std::endl;
        TR << "Error: Found residue position " << mp_conf.membrane() << " but looking for " << pose.conformation().chain_begin( nchains_final ) << std::endl;
        utility_exit_with_message("Integration test failed!");
    }

    // Check that every embedding residue in the embedding residue map actually maps to an EMB residue
    utility::vector1< std::pair< int, int > > embres_map = mp_conf.embres_map();
    for ( core::Size i = 1; i <= nchains_orig; ++i ) {

        core::Size resnum = embres_map[ i ].second;
        if ( pose.residue( resnum ).type().name().compare("EMB") != 0 ) {
            TR << "Error: Chain " << i << " embedding residue is incorrectly mapped!" << std::endl;
            TR << "Error: Looking for type EMB, but mapped embedding residue has type " << pose.residue( resnum ).type().name() << std::endl;
            utility_exit_with_message("Integration test failed!");
        }
    }

    // Check that the total residue in spanfile is consistent with the number of residues in the protein component of the pose
    // Also check number of chains is consistent with number of protein chains in the pose

    // Grab topologies from membrane conformation
    core::Size total_resnum_in_spanfiles = 0;
    utility::vector1< properties::SpanningTopology > topologies = mp_conf.spanning_topology();

    // compare number of chains
    if ( topologies.size() != nchains_orig ) {
        TR << "Error: Must specify a membrane spanning topology for every chain!" << std::endl;
        TR << "Error: Looking for " << nchains_orig << " topology objects but found " << topologies.size();
        utility_exit_with_message("Integration test failed!");
    }

    // COmpute total res
    for ( core::Size i = 1; i <= topologies.size(); ++i ) {
        core::Size nres = topologies[i].total_residue_in_span_file();
        total_resnum_in_spanfiles = total_resnum_in_spanfiles + nres;
    }

    // compare number of residues
    if ( total_resnum_in_spanfiles != pose.conformation().chain_end(nchains_orig) ) {
        TR << "Error: Total number of residues in spanning objects does not equal total number of residues in the pose" << std::endl;
        TR << "Error: Looking for " << pose.conformation().chain_end(nchains_orig) << " residues but found " << total_resnum_in_spanfiles << std::endl;
        utility_exit_with_message("Integration test failed!");
    }

    TR << "Test Passed!" << std::endl;

}

/// @brief      Test Membrane Fold Tree Topology
/// @details    Check that a valid memrbane fold tree has been constructed
void test_membrane_fold_tree( core::pose::Pose & pose ) {

    using namespace core::membrane;
    using namespace core::kinematics;

    TR << "Test: Membrane Fold Tree...";

    // Grab Membrane Conformation from the pose
    MembraneConformation & mp_conf = dynamic_cast< MembraneConformation & >( pose.conformation() );
    utility::vector1< std::pair< int, int > > embres_map = mp_conf.embres_map();

    // Derive the number of each type of jump (expected)
    core::Size chain_jumps = pose.conformation().num_chains()-2;
    core::Size membrane_jump = 1;
    core::Size embedding_jumps = pose.conformation().num_chains()-1;
    core::Size expected_total = chain_jumps + membrane_jump + embedding_jumps;

    // Check that all appropriate jump edges exist
    if ( pose.num_jump() != expected_total ) {
        TR << "Error: Number of jumps in the pose is not equal to the number of jumps expected for membrane foldtree" << std::endl;
        TR << "Error: Found " << pose.num_jump() << " but expected " << expected_total << std::endl;
        utility_exit_with_message("Integration test failed!");
    }

    // Check that the jump edges in the pose match up with that in conf
    FoldTree ft( pose.fold_tree() );
    for ( core::Size i = 1; i <= pose.conformation().num_chains()-1; ++i ) {
        if ( !ft.jump_exists( embres_map[ i ].first, embres_map[ i ].second ) ) {
            TR << "Error: Missing fold tree jump edge from position " << embres_map[ i ].first << " to " << embres_map[ i ].second << std::endl;
            TR << "Error: Or maybe your embedding residue map is constructed incorrectly?" << std::endl;
            utility_exit_with_message("Integraiton test failed!");
        }

    }

    // Check that the membrane jump exists
    if ( !ft.jump_exists( pose.conformation().chain_begin( pose.conformation().num_chains() ), 1 ) ) {
        TR << "Error: Missing membrane jump edge!" << std::endl;
        utility_exit_with_message("Integration test failed!");
    }

    // Checking that the membrane residue is the root of the foldtree
    if ( !ft.is_root( pose.conformation().chain_begin( pose.conformation().num_chains() ))) {
        TR << "Error: Membrane residue is not the root of the foldtree!" << std::endl;
        utility_exit_with_message("Integration test failed!");
    }

    // Checking that the membrane foldtree is a valid foldtree
    if ( !ft.check_fold_tree() ) {
        TR << "Error: Membrane fold tree is not a valid foldtree!" << std::endl;
        utility_exit_with_message("Integration test failed!");
    }

    TR << "Test Passed!" << std::endl;

}

/// @brief      Test Consistent Low Resolution MP Scoring
/// @details    Check that membrane scores do not change across iterations in lowres scoring
void test_membrane_scoring_lowres( core::pose::Pose & membrane_pose, core::pose::Pose & non_membrane_pose ) {

    using namespace core::scoring;
    using namespace protocols::simple_moves;

    TR << "Test: Low Resolution Membrane Scoring...";

    // Create a new Typeset Mover
    SwitchResidueTypeSetMoverOP typeset_mover = new SwitchResidueTypeSetMover( core::chemical::CENTROID );

    // Create Scoring Function
    core::scoring::ScoreFunctionOP lowres_membrane = core::scoring::ScoreFunctionFactory::create_score_function( "score_membrane" );

    // Switch residue typesets
    typeset_mover->apply(membrane_pose);
    typeset_mover->apply(non_membrane_pose);

    // Score Initial Object and Store pose energies
    lowres_membrane->score(membrane_pose);
    membrane_pose.energies().show(std::cout);

    lowres_membrane->score(non_membrane_pose);
    non_membrane_pose.energies().show(std::cout);

    // Calculate Score Energy Differences (sorry this is obnoxious)
    core::Real total_score = get_ddG( membrane_pose, non_membrane_pose, core::scoring::total_score );
    if ( total_score != 0 ) {
        TR << "Total score ddG != 0 from membrane lowres sfxn. Value is " << total_score << std::endl;
        utility_exit_with_message("Integration test failed!");
    }

    TR << "Test Passed!" << std::endl;
}

/// @brief      Test Consistent High Resolution MP Scoring
/// @details    Check that membrane scores do not change across iterations in high resolution scoring
void test_membrane_scoring_highres( core::pose::Pose & membrane_pose, core::pose::Pose & non_membrane_pose ) {

    using namespace core::scoring;

    TR << "Test: High Resolution Membrane Scoring...";

    // Create Scoring Function
    core::scoring::ScoreFunctionOP highres_membrane = core::scoring::ScoreFunctionFactory::create_score_function( "membrane_highres_Menv_smooth" );

    // Score Initial Object and Store pose energies
    highres_membrane->score(membrane_pose);
    highres_membrane->score(non_membrane_pose);

    // Calculate Score Energy Differences (sorry this is obnoxious)
    core::Real total_score = get_ddG( membrane_pose, non_membrane_pose, core::scoring::total_score );
    if ( total_score != 0 ) {
        TR << "Total score ddG != 0 from membrane highres sfxn. Value is " << total_score << std::endl;
        utility_exit_with_message("Integration test failed!");
    }

    TR << "Test Passed!" << std::endl;

}

/// @brief      Test Consistent Talaris Pose Scoring
/// @details    Check that membrane conformations do not change across iterations in standard
void test_talaris( core::pose::Pose & membrane_pose, core::pose::Pose & non_membrane_pose ) {

    using namespace core::scoring;

    TR << "Test: Talaris Scoring...";

    // Create Scoring Function
    core::scoring::ScoreFunctionOP talaris = core::scoring::ScoreFunctionFactory::create_score_function( "talaris2013" );

    // Score Initial Object and Store pose energies
    talaris->score(membrane_pose);
    talaris->score(non_membrane_pose);

    // Calculate Score Energy Differences (sorry this is obnoxious)
    core::Real total_score = get_ddG( membrane_pose, non_membrane_pose, core::scoring::total_score );
    if ( total_score != 0 ) {
        TR << "Total score ddG != 0 from talaris sfxn. Value is " << total_score << std::endl;
        utility_exit_with_message("Integration test failed!");
    }

    TR << "Test Passed!" << std::endl;
}

/// dump final pdb

/// @brief   Top Level Testing for Constructing Membrane Proteins
/// @details Runs applications with debug tracers contianing expected data. Will eventually
///          call the new object oriented code

/// @brief Main
int main( int argc, char* argv[] )
{

    try {

        using namespace protocols::membrane;
        using namespace protocols::jd2;

        // Initialize Options System, RG, and All Factory_Registrators
        devel::init(argc, argv);

        // Print App Info
        TR << "Integration Test: Membrane Protein Framework" << std::endl;
        TR << "Version: 1.0" << std::endl;
        TR << "Last Modified: 2/15/14" << std::endl;

        // Setup Pose //////////////////////////////////

        TR << "Setting up Membrane Protein from Initial Resources..." << std::endl;
        // Initialize Membrane Mover
        MembraneUnitTestMoverOP mp = new MembraneUnitTestMover();
        mp->register_options();
        mp->init_from_cmd();

        // Snd the mover through the job distributor
        JobDistributor::get_instance()->go(mp);

       // Grab Membrane Pose from the Membrane Mover
        core::pose::PoseOP membrane_pose = mp->get_membrane_pose();
        core::pose::PoseOP non_membrane_pose = mp->get_non_membrane_pose();

        TR << "Setup Complete!" << std::endl;

        // Run Benchamrking Tests //////////////////////////

        TR << "Running MP Framework Benchmark Tests" << std::endl;
       test_membrane_chain( *membrane_pose );
       test_membrane_conformation( *membrane_pose );
       test_membrane_fold_tree( *membrane_pose );
       test_talaris( *membrane_pose, *non_membrane_pose );

       // == tests currently not in use 2/15/14 ==
       //test_membrane_scoring_lowres( *membrane_pose, *non_membrane_pose );
       //test_membrane_scoring_highres( *membrane_pose, *non_membrane_pose );
       // test talaris lowres?

        TR << "Dumping the membrane pose to pdb" << std::endl;
        membrane_pose->dump_pdb("integration.pdb");

        TR << "All Tests Passed!" << std::endl;

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cout << "caught exception " << e.msg() << std::endl;
				return -1;
    }
}
