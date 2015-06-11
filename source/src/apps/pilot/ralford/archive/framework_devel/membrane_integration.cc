// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       apps/pilot/ralford/membrane_integration.cc
///
/// @brief      Top-Level Unit Test for the Membrane Protein Factory
/// @details    The purpose of this application is to test the membrane protein factory
///             initialization code including all external dependencies which cannot be tested
///             in JD2, Resource Manager, and the Pose cache. This can also serve as the integration
///             test for memrane protein initialization.
///
/// @note       This test is highly coupled to the xml file membrane.xml
/// @note       This test contains _no_ scoring integration!!!
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (1/2/14)

// App Headers
#include <devel/init.hh>

#include <protocols/membrane/MembraneMover.hh>

// Project Headers
#include <core/membrane/properties/SpanningTopology.hh>
#include <core/membrane/properties/MultiChainSpanningTopology.hh>

#include <core/membrane/util/definitions.hh>
#include <core/membrane/util/Exceptions.hh>

#include <core/membrane/kinematics/MembraneFoldTree.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/pose/membrane/util.hh>
#include <core/types.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

using basic::Error;
using basic::Warning;

using namespace core::membrane;

static thread_local basic::Tracer TR( "apps.pilot.ralford.membrane_integration" );

///// Individual Unit Tests ////////////

/// @brief      Test Contents of the Membrane Chain
/// @details    Check that the pose contains a membrane residue
///             and an embedding data residue for every original polymer chain
///             in the pose
void test_membrane_chain( core::pose::Pose & pose ) {

    // Calculate the number of chains in the original pose

    // Check the contents of the last chain
    // - first residue is a membrane residue
    // - the next n residues are embedding residues where n
    // is equal to the original n pose chains

}

/// @brief      Test the Contents of the Membrane Conformation
/// @details    Check that the position of the embedding residues is consistent
///             with teh numebring rules as well as the residue numbering for the membrane
///             residue. Check the contents of the embedding map is consistent with
///             membrane data. Check number of chains of topology and total res matches up
///             with the actual total residue.
void test_membrane_conformation( core::pose::Pose & pose ) {

    // check access info in membrane and embedding data getters

    // check consistency in the root

    // check consistency in the embedding data map

    // check that the total residue in span file (accross all span files) is equal to the total
    // number of residues in the pose minus the number of residues in the embedding
    // chain for consistency

    // check that the membrane conformation is valid

}

/// @brief      Test Membrane Fold Tree Topology
/// @details    Check that a valid memrbane fold tree has been constructed
void test_membrane_fold_tree( core::pose::Pose & pose ) {

    // Check that all appropriate jump edges exist

    // check the number of edges in the pose is correct

    // check the membrane residue is the root

    // check the fold tree is a valid fold tree

}

/// @brief      Test Consistent Low Resolution MP Scoring
/// @details    Check that membrane scores do not change across iterations in lowres scoring
void test_membrane_scoring_lowres( core::pose::Pose & pose ) {

    using namespace core::scoring;

    // calculate ddGs are 0

    //	pose.dump_scored_pdb( "repacked_once.pdb", *scorefxn );
//	TR << "Score after repacking once: " << pose.energies().total_energies()[ total_score ] << std::endl << std::endl;

}


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

        TR << "Starting to test the membrane mover!" << std::endl;

        // Initialize Membrane Mover
        MembraneMoverOP mp = new MembraneMover();
        JobDistributor::get_instance()->go(mp);

        // Grab Membrane Pose from the Membrane Mover
        core::pose::PoseOP pose = mp->get_pose();

        // Run Benchmarking tests
        test_membrane_chain( *pose );
        test_membrane_conformation( *pose );
        test_membrane_fold_tree( *pose );
        test_membrane_scoring_lowres( *pose ); // function I am going to write
        test_membrane_scoring_highres( *pose ): // function I am going to write
        test_talaris( *pose ); // function I am going to write

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cout << "caught exception " << e.msg() << std::endl;
				return -1;
    }
}
