// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/drug_design/SubstituentReplace.cxxtest.hh
/// @brief  test for SubstituentReplace Chemistry
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Project Headers
#include <core/types.hh>

#include <protocols/drug_design/SubstituentReplace.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/sdf/MolFileIOReader.hh>

// Utility Headers
#include <utility/file/FileName.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.drug_design.SubstituentReplace.cxxtest.hh");

// --------------- Test Class --------------- //

class SubstituentReplaceTests : public CxxTest::TestSuite {

private:
	core::chemical::ResidueTypeOP restype_;

public:

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	void test_substitute() {
		protocols::drug_design::SubstituentReplace replace;

		replace.template_database( "protocols/drug_design/pyrimidine_stub1.sdf");
		replace.substituents_database( "protocols/drug_design/ethoxypyrimidine.sdf");

		utility::vector1< core::chemical::MutableResidueTypeOP > restypes( core::chemical::sdf::convert_to_ResidueTypes("protocols/drug_design/carboxypyrimidine.sdf") );
		TS_ASSERT( restypes.size() > 0 );
		core::chemical::MutableResidueTypeOP restype( restypes[1] );

		core::chemical::NameVDMapping orig_map( *restype );

		TS_ASSERT_EQUALS( restype->natoms(), 13 );
		TS_ASSERT_EQUALS( restype->nheavyatoms(), 9 );

		//orig_map.show(TR);

		TS_ASSERT( restype->has("O1") );
		TS_ASSERT( restype->has("O2") );
		TS_ASSERT( restype->has("C5") );
		TS_ASSERT( restype->has("N1") );
		TS_ASSERT( restype->has("N2") );

		replace.apply(*restype);
		TS_ASSERT_EQUALS( replace.get_last_status(), core::chemical::modifications::SUCCESS );

		TS_ASSERT_EQUALS( restype->natoms(), 17 );
		TS_ASSERT_EQUALS( restype->nheavyatoms(), 10 );

		// Keeps the nitrogens
		TS_ASSERT( restype->has("N1") );
		TS_ASSERT( restype->has("N2") );
		TS_ASSERT( restype->has("O1") ); // Has the one oxygen
		TS_ASSERT( ! restype->has("O2") ); // Has the one oxygen

		TS_ASSERT( restype->has("F1") ); // Adds the flurine
		TS_ASSERT( ! restype->has("CL1") ); // Does not add the chlorine

		core::chemical::VDVDMapping map( replace.get_mapping() );

		TS_ASSERT_EQUALS( map[ orig_map[ " O2 " ] ], map.invalid_entry() );

		TS_ASSERT_EQUALS( restype->atom_name( map[ orig_map[ " N1 " ] ] ), " N1 " );
		TS_ASSERT_EQUALS( restype->atom_name( map[ orig_map[ " N2 " ] ] ), " N2 " );

		// Shouldn't map anything on the flourine or oxygen
		TS_ASSERT_EQUALS( map.reverse_lookup( restype->atom_vertex("F1") ), map.invalid_key() );
		TS_ASSERT_EQUALS( map.reverse_lookup( restype->atom_vertex("O1") ), map.invalid_key() );
	}

	// Now go the other direction
	void test_substitute_rev() {
		protocols::drug_design::SubstituentReplace replace;

		replace.template_database( "protocols/drug_design/pyrimidine_stub1.sdf");
		replace.substituents_database( "protocols/drug_design/carboxypyrimidine.sdf");

		utility::vector1< core::chemical::MutableResidueTypeOP > restypes( core::chemical::sdf::convert_to_ResidueTypes("protocols/drug_design/ethoxypyrimidine.sdf") );
		TS_ASSERT( restypes.size() > 0 );
		core::chemical::MutableResidueTypeOP restype( restypes[1] );

		core::chemical::NameVDMapping orig_map( *restype );

		TS_ASSERT_EQUALS( restype->natoms(), 17 );
		TS_ASSERT_EQUALS( restype->nheavyatoms(), 11 );

		//orig_map.show(TR);

		TS_ASSERT( restype->has("O1") );
		TS_ASSERT( restype->has("N1") );
		TS_ASSERT( restype->has("N2") );
		TS_ASSERT( restype->has("F1") );
		TS_ASSERT( restype->has("CL1") );

		replace.apply(*restype);
		TS_ASSERT_EQUALS( replace.get_last_status(), core::chemical::modifications::SUCCESS );

		TS_ASSERT_EQUALS( restype->natoms(), 13 );
		TS_ASSERT_EQUALS( restype->nheavyatoms(), 10 );

		// Keeps the nitrogens
		TS_ASSERT( restype->has("N1") );
		TS_ASSERT( restype->has("N2") );
		TS_ASSERT( restype->has("O1") ); // Has the two oxygens
		TS_ASSERT( restype->has("O2") );
		TS_ASSERT( restype->has("CL1") ); // keeps the chlorine
		TS_ASSERT( ! restype->has("F1") ); // loses the flurine

		core::chemical::VDVDMapping map( replace.get_mapping() );

		TS_ASSERT_EQUALS( map[ orig_map[ " F1 " ] ], map.invalid_entry() );
		TS_ASSERT_EQUALS( map[ orig_map[ " O1 " ] ], map.invalid_entry() );

		TS_ASSERT_EQUALS( restype->atom_name( map[ orig_map[ " N1 " ] ] ), " N1 " );
		TS_ASSERT_EQUALS( restype->atom_name( map[ orig_map[ " N2 " ] ] ), " N2 " );
		TS_ASSERT_EQUALS( restype->atom_name( map[ orig_map[ " CL1" ] ] ), " CL1" );

		// Shouldn't map anything on the oxygens
		TS_ASSERT_EQUALS( map.reverse_lookup( restype->atom_vertex("O1") ), map.invalid_key() );
		TS_ASSERT_EQUALS( map.reverse_lookup( restype->atom_vertex("O2") ), map.invalid_key() );
	}

};
