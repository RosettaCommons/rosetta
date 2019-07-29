// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/ResidueTypeSetTests.cxxtest.hh
/// @brief unit test for ResidueTypeFinder
/// @author Rocco Moretti (rmorettiase@gmail.com),

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>

// Project Headers
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>

// Platform Headers
#include <basic/Tracer.hh>
#include <basic/database/open.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ Headers
#include <string>
#include <ostream>

static basic::Tracer TR("core.chemical.ResidueTypeFinderTests.cxxtest");

class ResidueTypeFinderTests : public CxxTest::TestSuite {

	core::chemical::ResidueTypeSetCOP rts_;

public:

	void setUp() {
		using namespace core::chemical;
		core_init();
		rts_ = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	}

	void tearDown() {}

	void test_get_by_aa() {
		using namespace core::chemical;

		ResidueTypeFinder rtf( *rts_ );
		rtf.aa( aa_ala );
		ResidueTypeCOPs residues( rtf.get_possible_base_residue_types() );
		TS_ASSERT_EQUALS( residues.size(), 1 );
	}

	void test_get_by_std_name3() {
		using namespace core::chemical;

		ResidueTypeFinder rtf( *rts_ );
		rtf.name3( "ALA" );
		ResidueTypeCOPs residues( rtf.get_possible_base_residue_types() );
		TS_ASSERT_EQUALS( residues.size(), 1 );
		TR << "Redoing " << std::endl;
		// Check to make sure we don't accidentally pollute with CCD on second time around.
		ResidueTypeFinder rtf2( *rts_ );
		rtf2.name3( "ALA" );
		ResidueTypeCOPs residues2( rtf2.get_possible_base_residue_types() );
		TS_ASSERT_EQUALS( residues2.size(), 1 );
		TR << "Redoing done " << std::endl;
	}

	void test_get_by_odd_name3() {
		using namespace core::chemical;

		ResidueTypeFinder rtf( *rts_ );
		{
			rtf.name3( "AZ&" );
			ResidueTypeCOPs residues( rtf.get_possible_base_residue_types() );
			TS_ASSERT_EQUALS( residues.size(), 0 );
		}

		{
			rtf.name3( "AZT" );
			ResidueTypeCOPs residues( rtf.get_possible_base_residue_types() );
			TS_ASSERT_EQUALS( residues.size(), 1 ); // now in components.cif
			TS_ASSERT_EQUALS( residues[1]->name(), "pdb_AZT" );
		}
	}

	void test_base_overlapping_name3() {
		using namespace core::chemical;

		ResidueTypeFinder rtf( *rts_ );
		rtf.name3( "C12" );
		{
			ResidueTypeCOPs residues( rtf.get_possible_base_residue_types() );
			TS_ASSERT_EQUALS( residues.size(), 3 ); // database C12, D-patched version and CCD version
		}
		// Check that we don't double-add CCD to base type list
		{
			ResidueTypeCOPs residues( rtf.get_possible_base_residue_types() );
			TS_ASSERT_EQUALS( residues.size(), 3 ); // database C12, D-patched version and CCD version
		}
	}

	void test_base_overlapping_name3_no_CCD() {
		using namespace core::chemical;

		ResidueTypeFinder rtf( *rts_ );
		rtf.set_no_CCD_on_name3_match( true );
		rtf.name3( "C12" );
		ResidueTypeCOPs residues( rtf.get_possible_base_residue_types(/*include_unpatchable*/ true, /*apply_all_filters*/ true) );
		// need apply_all_filters true to apply the no_CCD filtering.
		// While we won't necessarily load the CCD version, if it's already been loaded it would show up without apply_all_filters
		TS_ASSERT_EQUALS( residues.size(), 2 ); // database C12, D-patched version and no CCD version
	}


	void test_get_by_all_possible() {
		using namespace core::chemical;

		ResidueTypeFinder rtf( *rts_ );
		rtf.name3( "C12" );
		{
			ResidueTypeCOPs residues( rtf.get_all_possible_residue_types() );
			TS_ASSERT_EQUALS( residues.size(), 3 ); // database C12, D-patched version and CCD version
		}
		// Check we don't double-add CCD version
		{
			ResidueTypeCOPs residues( rtf.get_all_possible_residue_types() );
			TS_ASSERT_EQUALS( residues.size(), 3 ); // database C12, D-patched version and CCD version
		}
	}

	void test_get_by_all_possible_no_CCD() {
		using namespace core::chemical;

		ResidueTypeFinder rtf( *rts_ );
		rtf.set_no_CCD_on_name3_match( true );
		rtf.name3( "C12" );
		ResidueTypeCOPs residues( rtf.get_all_possible_residue_types() );
		TS_ASSERT_EQUALS( residues.size(), 2 ); // database C12, D-patched version and no CCD version
	}

	void test_get_best_atom_names() {
		using namespace core::chemical;

		ResidueTypeFinder rtf( *rts_ );
		rtf.name3( "C12" );
		ResidueTypeCOP res1 = rtf.get_best_match_residue_type_for_atom_names( { "N", "CA", "CB", "CE" } ); // Names from database version
		TS_ASSERT_EQUALS( res1->name(), "C12" );
		ResidueTypeCOP res2 = rtf.get_best_match_residue_type_for_atom_names( { "C1", "N2", "CG1", "OG1", "C2" } ); // Names from CCD version
		TS_ASSERT_EQUALS( res2->name(), "pdb_C12" );
	}

	void test_get_best_atom_names_no_CCD() {
		using namespace core::chemical;

		ResidueTypeFinder rtf( *rts_ );
		rtf.set_no_CCD_on_name3_match( true );
		rtf.name3( "C12" );
		ResidueTypeCOP res1 = rtf.get_best_match_residue_type_for_atom_names( { "N", "CA", "CB", "CE" } ); // Names from database version
		TS_ASSERT_EQUALS( res1->name(), "C12" );
		ResidueTypeCOP res2 = rtf.get_best_match_residue_type_for_atom_names( { "C1", "N2", "CG1", "OG1", "C2" } ); // Names from CCD version
		TS_ASSERT_EQUALS( res2->name(), "C12" ); // Matches zero names, but that's what the set_no_CCD_on_name3_match(true) implies ...
	}

};
