// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/chemical/CacheableResidueTypeSetsTests.cxxtest.hh
/// @brief  Test CacheableResidueTypeSets
/// @author Rocco Moretti (rmorettiase@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

// Project Headers
#include <core/chemical/CacheableResidueTypeSets.hh>

// Core Headers
#include <core/chemical/PoseResidueTypeSet.hh>

#include <basic/Tracer.hh>

#include <core/chemical/ChemicalManager.hh> // AUTO IWYU For ChemicalManager
#include <core/init_util.hh> // AUTO IWYU For core_init

static basic::Tracer TR("CacheableResidueTypeSetsTests.cxxtest");


class CacheableResidueTypeSetsTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	void test_setting() {

		using namespace core::chemical;

		core::chemical::CacheableResidueTypeSetsOP cache( new core::chemical::CacheableResidueTypeSets() );

		// We start off without the RTS, and we get a nullptr if we try to access them.
		TS_ASSERT( ! cache->has_res_type_set( FULL_ATOM_t ) );
		TS_ASSERT( ! cache->has_res_type_set( CENTROID_t ) );
		TS_ASSERT_EQUALS( cache->get_res_type_set( FULL_ATOM_t ), nullptr );
		TS_ASSERT_EQUALS( cache->get_res_type_set( CENTROID_t ), nullptr );

		PoseResidueTypeSetOP fa( new PoseResidueTypeSet( ChemicalManager::get_instance()->residue_type_set( FULL_ATOM_t ) ) );
		PoseResidueTypeSetOP cen( new PoseResidueTypeSet( ChemicalManager::get_instance()->residue_type_set( CENTROID_t ) ) );

		cache->set_res_type_set( cen );
		TS_ASSERT( ! cache->has_res_type_set( FULL_ATOM_t ) );
		TS_ASSERT( cache->has_res_type_set( CENTROID_t ) );

		cache->set_res_type_set( fa, FULL_ATOM_t );
		TS_ASSERT( cache->has_res_type_set( FULL_ATOM_t ) );
		TS_ASSERT( cache->has_res_type_set( CENTROID_t ) );

		// We get the right mode back when we call the respective:
		TS_ASSERT_EQUALS( cache->get_res_type_set( FULL_ATOM_t )->mode(), FULL_ATOM_t );
		TS_ASSERT_EQUALS( cache->get_res_type_set( CENTROID_t )->mode(), CENTROID_t );

		// Test that we can set them to nullptrs
		cache->set_res_type_set( nullptr, CENTROID_t );
		TS_ASSERT( cache->has_res_type_set( FULL_ATOM_t ) );
		TS_ASSERT( ! cache->has_res_type_set( CENTROID_t ) );

		// We can force a RTS into the wrong slot, though we'll get a message about it.
		// (This isn't particularly useful, so don't let this unit test keep you from
		// changing this behavior.)
		cache->set_res_type_set( cen, FULL_ATOM_t );
		TS_ASSERT( cache->has_res_type_set( FULL_ATOM_t ) );
		TS_ASSERT( ! cache->has_res_type_set( CENTROID_t ) );

		TS_ASSERT_EQUALS( cache->get_res_type_set( FULL_ATOM_t )->mode(), CENTROID_t ); // Stupid, but that's the way it is.

	}

};



