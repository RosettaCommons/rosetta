// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @file membrane_pose_test.cc (quick debug)
///
/// @brief Quick debugging of constructing a membrane pose (enables gdb stuff)
/// @details
///
/// @author Rebecca Alford (rfalford12@gmail.com)

// Package Headers
#include <core/pose/Pose.hh>

#include <devel/init.hh>
#include <core/types.hh>

// Project Headers
#include <core/membrane/util/definitions.hh>
#include <core/membrane/util/Exceptions.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>

#include <basic/Tracer.hh>

// Tested Classes
#include <core/membrane/geometry/MembraneResidueFactory.hh>
#include <core/membrane/geometry/EmbeddingFactory.hh>
#include <core/membrane/io/EmbedDefIO.hh>
#include <core/membrane/io/SpanFileIO.hh>

#include <core/membrane/kinematics/MembraneFoldTree.hh>

// Package headers
#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/string_util.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

using basic::Error;
using basic::Warning;

static basic::Tracer TR( "apps.pilot.ralford.membrane_pose_test" );

/// @brief   Top Level Testing for Constructing Membrane Proteins
/// @details Runs applications with debug tracers contianing expected data. Will eventually
///          call the new object oriented code

/// @brief Main
int main( int argc, char* argv[] )
{
    
    using namespace core::chemical;
    using namespace core::conformation;
    
    using namespace basic::options;
    using namespace core::membrane::geometry;
    using namespace core::membrane::kinematics;
    using namespace core::pose;
    using namespace core::import_pose;
    using namespace core::conformation;
    using namespace core::membrane;
    
    // Initialize
    devel::init( argc, argv );
    
    core::pose::PoseOP pose = new Pose();
    pose_from_pdb( *pose, "test/core/membrane/io/1afo_A.pdb");
    
    // Membrane Residue Factory
    MembraneResidueFactory mrf;
    
    // Create and add Membrane residue with defaults
    core::Vector center(0, 0, 0);
    core::Vector normal(0, 0, 1);
    core::Real depth = 30.0;
    
    mrf.add_membrane_residue(center, normal, depth, *pose, true);
    mrf.add_embedding_residue(center, normal, depth, *pose, 1, true);
    
    TR << "Showing foldtree after pose construction" << std::endl;
    pose->fold_tree().show(std::cout);\
    
    TR << "Checking what I know - is this a fold tree?" << std::endl;
    if ( pose->fold_tree().check_fold_tree() ) {
        TR << "YES" << std::endl;
    }
    
    TR << "Going to try to construct a foldtree from this one" << std::endl;
    std::map< int, int > chain_map1;
    chain_map1.insert( std::pair< int, int >( 1, 42 ) );
    
    MembraneFoldTreeOP ft = new MembraneFoldTree( pose->fold_tree(), chain_map1, 41 );
    
    TR << "Showing new foldtree" << std::endl;
    ft->show(std::cout);
    
    TR << "The big question, is it a fold tree??" << std::endl;
    if ( ft->check_fold_tree() ) {
        TR << "YES" << std::endl;
    }
    
    return 0;
}