// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    test/core/chemical/util.cxxtest.hh
/// @brief   Test suite for chemistry-related utility functions
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <core/chemical/util.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>

#include <basic/Tracer.hh>


static THREAD_LOCAL basic::Tracer TR( "core.chemical.util.cxxtest.hh" );

class ChemistryNamespaceUtilityFunctionTests : public CxxTest::TestSuite {
public:  // Standard methods //////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		using namespace core::chemical;

		core_init_with_additional_options( "-include_sugars" );

		ChemicalManager * manager( ChemicalManager::get_instance() );
		type_set_ = manager->residue_type_set( FA_STANDARD );
	}

	// Destruction
	void tearDown()
	{}


public:  // Tests /////////////////////////////////////////////////////////////
	// Confirm that Rosetta can figure out if a main-chain torsion is also a ring torsion.
	void test_is_mainchain_torsion_also_ring_torsion()
	{
		using namespace core;
		using namespace core::chemical;

		ResidueType const & alanine( type_set_->name_map( "ALA" ) );
		ResidueType const & to4_glucose( type_set_->name_map( "->4)-alpha-D-Glcp" ) );
		ResidueType const & to6_glucose( type_set_->name_map( "->6)-alpha-D-Glcp" ) );

		// Alanine is acyclic, so none of the main-chain torsions should be ring torsions.
		Size const n_Ala_mainchain_torsions( alanine.mainchain_atoms().size() );
		for ( core::uint i( 1 ); i <= n_Ala_mainchain_torsions; ++ i ) {
			TS_ASSERT( ! is_mainchain_torsion_also_ring_torsion( alanine, i ) );
		}

		// ->4)-alpha-D-Glcp has 5 main-chain torsions.
		TS_ASSERT( is_mainchain_torsion_also_ring_torsion( to4_glucose, 1 ) );  // also nu1
		TS_ASSERT( is_mainchain_torsion_also_ring_torsion( to4_glucose, 2 ) );  // also nu2
		TS_ASSERT( is_mainchain_torsion_also_ring_torsion( to4_glucose, 3 ) );  // also nu3
		TS_ASSERT( ! is_mainchain_torsion_also_ring_torsion( to4_glucose, 4 ) );  // psi(n+1)
		TS_ASSERT( ! is_mainchain_torsion_also_ring_torsion( to4_glucose, 5 ) );  // phi(n+1)

		// ->6)-alpha-D-Glcp has 7 main-chain torsions.
		TS_ASSERT( is_mainchain_torsion_also_ring_torsion( to6_glucose, 1 ) );  // also nu1
		TS_ASSERT( is_mainchain_torsion_also_ring_torsion( to6_glucose, 2 ) );  // also nu2
		TS_ASSERT( is_mainchain_torsion_also_ring_torsion( to6_glucose, 3 ) );  // also nu3
		TS_ASSERT( is_mainchain_torsion_also_ring_torsion( to6_glucose, 4 ) );  // also nu4
		TS_ASSERT( ! is_mainchain_torsion_also_ring_torsion( to6_glucose, 5 ) );  // psi(n+1)
		TS_ASSERT( ! is_mainchain_torsion_also_ring_torsion( to6_glucose, 6 ) );  // phi(n+1)
		TS_ASSERT( ! is_mainchain_torsion_also_ring_torsion( to6_glucose, 7 ) );  // omega(n+1)
	}

private:  // Private data /////////////////////////////////////////////////////
	core::chemical::ResidueTypeSetCOP type_set_;

};  // class ChemistryNamespaceUtilityFunctionTests
