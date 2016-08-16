// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file membrane_database_checks.cc
///
/// @brief   Quick Checking of membrane residue and atom types
/// @details This will eventually be converted into an appropriate db io unit test
///
/// @author Rebecca Alford (rfalford12@gmail.com)

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <core/types.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/AA.hh>

#include <core/chemical/AtomType.hh>

// Utility Headers
#include <utility/vector1.hh> 
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib> 
#include <string> 

using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.ralford.mp_database_check" );

/// @brief Load in a Centroid Pose and Check rsd/atom types
void
check_centroid() {

    using namespace core::chemical;
    using namespace core::conformation;

    TR << "======Testing membrane protein centroid database components======" << std::endl;

    // Residue type information
    ResidueTypeSetCAP const & residue_set(
                                          core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID )
                                          );

    // Create a membrane and embedding residue
    ResidueTypeCOPs const & emb_list( residue_set->name3_map("EMB") );
    ResidueTypeCOPs const & mem_list( residue_set->name3_map("MEM") );

    ResidueType const & membrane( *mem_list[1] );
    ResidueType const & embedding( *emb_list[1] );

    ResidueOP mem( ResidueFactory::create_residue(membrane) );
    ResidueOP emb( ResidueFactory::create_residue(embedding) );

    // Get AA base types
    AA mem_aa = mem->aa();
    AA emb_aa = emb->aa();
    
    // Print Residue properties
    TR << "Membrane Residue properties" << std::endl;
    TR << "Name: " << mem->name() << std::endl;
    TR << "AA Base Type: " << name_from_aa( mem_aa ) << std::endl;
    TR << "Atom Type: " << mem->type().atom_type( 1 ).atom_type_name() << std::endl;
    TR << "Atom Type: " << mem->type().atom_type( 2 ).atom_type_name() << std::endl;
    TR << "Atom Type: " << mem->type().atom_type( 3 ).atom_type_name() << std::endl;

    TR << "Embedding Residue Properties" << std::endl;
    TR << "Name: " << emb->name() << std::endl;
    TR << "AA Base Type: " << name_from_aa( emb_aa ) << std::endl;
    TR << "Atom Type: " << emb->type().atom_type( 1 ).atom_type_name() << std::endl;
    TR << "Atom Type: " << emb->type().atom_type( 2 ).atom_type_name() << std::endl;
    TR << "Atom Type: " << emb->type().atom_type( 3 ).atom_type_name() << std::endl;
    TR << "Centroid DB Test complete!" << std::endl;

}

/// @brief Load in a Fullatom Pose and Check rsd/atom types
void
check_fullatom() {

	using namespace core::chemical;
    using namespace core::conformation;

    TR << "======Testing membrane protein fullatom database components======" << std::endl;

        // Residue type information
    ResidueTypeSetCAP const & residue_set(
                                          core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )
                                          );

    // Create a membrane and embedding residue
    ResidueTypeCOPs const & emb_list( residue_set->name3_map("EMB") );
    ResidueTypeCOPs const & mem_list( residue_set->name3_map("MEM") );

    ResidueType const & membrane( *mem_list[1] );
    ResidueType const & embedding( *emb_list[1] );

    ResidueOP mem( ResidueFactory::create_residue(membrane) );
    ResidueOP emb( ResidueFactory::create_residue(embedding) );

    // Get AA base types
    AA mem_aa = mem->aa();
    AA emb_aa = emb->aa();
    
    // Print Residue properties
    TR << "Membrane Residue properties" << std::endl;
    TR << "Name: " << mem->name() << std::endl;
    TR << "AA Base Type: " << name_from_aa( mem_aa ) << std::endl;
    TR << "Atom Type: " << mem->type().atom_type( 1 ).atom_type_name() << std::endl;
    TR << "Atom Type: " << mem->type().atom_type( 2 ).atom_type_name() << std::endl;
    TR << "Atom Type: " << mem->type().atom_type( 3 ).atom_type_name() << std::endl;

    TR << "Embedding Residue Properties" << std::endl;
    TR << "Name: " << emb->name() << std::endl;
    TR << "AA Base Type: " << name_from_aa( emb_aa ) << std::endl;
    TR << "Atom Type: " << emb->type().atom_type( 1 ).atom_type_name() << std::endl;
    TR << "Atom Type: " << emb->type().atom_type( 2 ).atom_type_name() << std::endl;
    TR << "Atom Type: " << emb->type().atom_type( 3 ).atom_type_name() << std::endl;

    TR << "Fullatom DB Test complete!" << std::endl;
    
}

/// @brief Main
int main( int argc, char* argv[] )
{
    
    using namespace core::chemical;
    using namespace core::conformation;
    
    // Initialize
    devel::init( argc, argv );
    
    TR << "Checking membrane protein database properties" << std::endl;

    // Run for both centroid and fullatom
    check_centroid();
    check_fullatom();

    TR << "Test complete!" << std::endl;

}

