// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/MPDatabaseIO.cxxtest.hh
///
/// @brief 		Test Suite for Membrane Protein Derived Atom and Residue Types
/// @details    CxxTest suite
///
/// @author 	Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified: 1/26/14

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

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

class MPDatabaseIOTest : public CxxTest::TestSuite {
    
public: // test functions
    
    /// Test Setup Functions ////////
    
    /// @brief Setup Test
    void setUp() {
        
        // Initialize
        core_init();
    }
    
    /// @brief Standard Tear Down
    void tearDown() {}


    /// @brief Load in a Centroid Pose and Check rsd/atom types
    void test_centroid_db() {
    
        using namespace core::chemical;
        using namespace core::conformation;
        
        TS_TRACE("Testing membrane protein centroid database attributes...");
    
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
        
        // String type checking for membrane residues
        TS_ASSERT_EQUALS( mem->name().compare("MEM"), 0);
        TS_ASSERT_EQUALS( name_from_aa( mem_aa ).compare("MPR"), 0);
        TS_ASSERT_EQUALS( mem->type().atom_type( 1 ).atom_type_name().compare("MPnm"), 0 );
        TS_ASSERT_EQUALS( mem->type().atom_type( 2 ).atom_type_name().compare("MPct"), 0 );
        TS_ASSERT_EQUALS( mem->type().atom_type( 3 ).atom_type_name().compare("MPtk"), 0 );
        
        // String type checking for embedding residues
        TS_ASSERT_EQUALS( emb->name().compare("EMB"), 0);
        TS_ASSERT_EQUALS( name_from_aa( emb_aa ).compare("MPR"), 0);
        TS_ASSERT_EQUALS( emb->type().atom_type( 1 ).atom_type_name().compare("MPnm"), 0 );
        TS_ASSERT_EQUALS( emb->type().atom_type( 2 ).atom_type_name().compare("MPct"), 0 );
        TS_ASSERT_EQUALS( emb->type().atom_type( 3 ).atom_type_name().compare("MPdp"), 0 );
    
    }

    /// @brief Load in a Fullatom Pose and Check rsd/atom types
    void test_fullatom_db() {
    
        using namespace core::chemical;
        using namespace core::conformation;
        
        TS_TRACE("Testing membrane protein fullatom (fa_standard) database attributes...");
        
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
        
        // String type checking for membrane residues
        TS_ASSERT_EQUALS( mem->name().compare("MEM"), 0);
        TS_ASSERT_EQUALS( name_from_aa( mem_aa ).compare("MPR"), 0);
        TS_ASSERT_EQUALS( mem->type().atom_type( 1 ).atom_type_name().compare("MPnm"), 0 );
        TS_ASSERT_EQUALS( mem->type().atom_type( 2 ).atom_type_name().compare("MPct"), 0 );
        TS_ASSERT_EQUALS( mem->type().atom_type( 3 ).atom_type_name().compare("MPtk"), 0 );
        
        // String type checking for embedding residues
        TS_ASSERT_EQUALS( emb->name().compare("EMB"), 0);
        TS_ASSERT_EQUALS( name_from_aa( emb_aa ).compare("MPR"), 0);
        TS_ASSERT_EQUALS( emb->type().atom_type( 1 ).atom_type_name().compare("MPnm"), 0 );
        TS_ASSERT_EQUALS( emb->type().atom_type( 2 ).atom_type_name().compare("MPct"), 0 );
        TS_ASSERT_EQUALS( emb->type().atom_type( 3 ).atom_type_name().compare("MPdp"), 0 );
    
    }
}; // test suite - MPDatabase IO

