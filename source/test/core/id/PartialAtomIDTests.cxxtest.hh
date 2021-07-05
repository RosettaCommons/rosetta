// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/id/PartialAtomIDTests.cxxtest.hh
/// @brief  Tests for the partial atom ID class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

// Project Headers
#include <core/id/PartialAtomID.hh>

// Core Headers

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <sstream>

static basic::Tracer TR("PartialAtomIDTests");


class PartialAtomIDTests : public CxxTest::TestSuite {
	//Define Variables

public:
	typedef core::id::PartialAtomID PartialAtomID;

	void setUp(){}

	void tearDown(){}

	void test_partial_atom_id_default_constructor() {
		PartialAtomID id;
		TS_ASSERT_EQUALS( id.rsd(), 0 );
		TS_ASSERT_EQUALS( id.atomno(), 0 );
		TS_ASSERT_EQUALS( id.valid(), false );
		TS_ASSERT_EQUALS( id.complete(), false );
		TS_ASSERT_EQUALS( id.partial(), false );
	}

	void test_construct_fully_resolved_partial_atom_id() {
		PartialAtomID id(5, 22);
		TS_ASSERT_EQUALS( id.atomno(), 5 );
		TS_ASSERT_EQUALS( id.rsd(), 22 );
		TS_ASSERT_EQUALS( id.valid(), true );
		TS_ASSERT_EQUALS( id.complete(), true );
		TS_ASSERT_EQUALS( id.partial(), false );
	}

	void test_construct_partially_resolved_partial_atom_id() {
		// the upper-connect atom on residue 22, assuming it has a lower- and upper connect
		PartialAtomID id(2, 22, 0);
		TS_ASSERT_EQUALS( id.resconnid(), 2 );
		TS_ASSERT_EQUALS( id.rsd(), 22 );
		TS_ASSERT_EQUALS( id.bonds_from_resconn(), 0 );
		TS_ASSERT_EQUALS( id.valid(), true );
		TS_ASSERT_EQUALS( id.complete(), false );
		TS_ASSERT_EQUALS( id.partial(), true );
	}

	void test_partial_atom_id_set_partial() {
		PartialAtomID id( 5, 22 ); // fully resolved
		id.set_partial( 2, 0, 22 );
		TS_ASSERT_EQUALS( id.resconnid(), 2 );
		TS_ASSERT_EQUALS( id.rsd(), 22 );
		TS_ASSERT_EQUALS( id.bonds_from_resconn(), 0 );
		TS_ASSERT_EQUALS( id.valid(), true );
		TS_ASSERT_EQUALS( id.complete(), false );
		TS_ASSERT_EQUALS( id.partial(), true );
	}

	void test_partial_atom_id_set_complete() {
		PartialAtomID id( 2, 22, 0 );
		id.set_complete(5, 22);
		TS_ASSERT_EQUALS( id.atomno(), 5 );
		TS_ASSERT_EQUALS( id.rsd(), 22 );
		TS_ASSERT_EQUALS( id.valid(), true );
		TS_ASSERT_EQUALS( id.complete(), true );
		TS_ASSERT_EQUALS( id.partial(), false );
	}

	void test_partial_atom_id_invalid_from_full_ctor() {
		PartialAtomID id( 0, 22);
		TS_ASSERT_EQUALS( id.atomno(), 0 );
		TS_ASSERT_EQUALS( id.rsd(), 22 );
		TS_ASSERT_EQUALS( id.valid(), false );
		TS_ASSERT_EQUALS( id.complete(), false );
		TS_ASSERT_EQUALS( id.partial(), false );
	}

	void test_partial_atom_id_invalid_from_partial_ctor() {
		PartialAtomID id( 0, 22, 0);
		TS_ASSERT_EQUALS( id.resconnid(), 0 );
		TS_ASSERT_EQUALS( id.bonds_from_resconn(), 0 );
		TS_ASSERT_EQUALS( id.rsd(), 22 );
		TS_ASSERT_EQUALS( id.valid(), false );
		TS_ASSERT_EQUALS( id.complete(), false );
		TS_ASSERT_EQUALS( id.partial(), false );
	}

	void test_partial_atom_id_stream_output() {
		PartialAtomID id1( 5, 22 );
		std::ostringstream oss1;
		oss1 << id1;
		TS_ASSERT_EQUALS( oss1.str(), "atomno= 5 resconn_id= 0 bonds_from_resconn= 0 rsd= 22" );

		PartialAtomID id2( 2, 22, 0 );
		std::ostringstream oss2;
		oss2 << id2;
		TS_ASSERT_EQUALS( oss2.str(), "atomno= 0 resconn_id= 2 bonds_from_resconn= 0 rsd= 22" );
	}

	void test_partial_atom_id_ordering() {
		PartialAtomID id1( 1, 22 );
		PartialAtomID id2( 2, 22 );
		PartialAtomID id3( 0, 23, 0);
		PartialAtomID id4( 0, 23, 1);
		PartialAtomID id5( 1, 23, 0);
		PartialAtomID id6( 1, 23, 2);
		PartialAtomID id7( 1, 23 );
		PartialAtomID id8( 2, 23 );

		std::vector< PartialAtomID > ids( {id1, id2, id3, id4, id5, id6, id7, id8} );

		for ( core::Size ii = 0; ii < ids.size(); ++ii ) {
			for ( core::Size jj = ii+1; jj < ids.size(); ++jj ) {
				TS_ASSERT_LESS_THAN( ids[ii], ids[jj] );
			}
		}
	}

};
