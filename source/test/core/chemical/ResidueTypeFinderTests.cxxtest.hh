// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
	}

	void test_get_by_odd_name3() {
		using namespace core::chemical;

		ResidueTypeFinder rtf( *rts_ );
		rtf.name3( "AZT" );
		ResidueTypeCOPs residues( rtf.get_possible_base_residue_types() );
		TS_ASSERT_EQUALS( residues.size(), 0 );
	}

	void test_get_by_overlapping_name3() {
		using namespace core::chemical;

		ResidueTypeFinder rtf( *rts_ );
		rtf.name3( "C12" );
		ResidueTypeCOPs residues( rtf.get_possible_base_residue_types() );
		TS_ASSERT_EQUALS( residues.size(), 2 ); // l-aa and d-aa
	}
};

class ResidueTypeFinderComponentsTests : public CxxTest::TestSuite {

	core::chemical::ResidueTypeSetCOP rts_;

public:

	void setUp() {
		using namespace core::chemical;
		core_init_with_additional_options("-load_PDB_components -PDB_components_file core/chemical/mmCIF/components_trimmed.cif");
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
		TS_ASSERT_EQUALS( residues.size(), 1 ); // This is likely to fail now.
		std::cout << "Redoing " << std::endl;
		ResidueTypeFinder rtf2( *rts_ );
		rtf2.name3( "ALA" );
		ResidueTypeCOPs residues2( rtf2.get_possible_base_residue_types() );
		TS_ASSERT_EQUALS( residues2.size(), 1 ); // Should be quick, as ALA_pdb is prohibited
		std::cout << "Redoing done " << std::endl;
	}

	void test_get_by_odd_name3() {
		using namespace core::chemical;

		ResidueTypeFinder rtf( *rts_ );
		rtf.name3( "AZT" );
		ResidueTypeCOPs residues( rtf.get_possible_base_residue_types() );
		TS_ASSERT_EQUALS( residues.size(), 1 );
	}

	void test_get_by_overlapping_name3() {
		using namespace core::chemical;

		ResidueTypeFinder rtf( *rts_ );
		rtf.name3( "C12" );
		ResidueTypeCOPs residues( rtf.get_possible_base_residue_types() );
		TS_ASSERT_EQUALS( residues.size(), 3 );
	}
};
