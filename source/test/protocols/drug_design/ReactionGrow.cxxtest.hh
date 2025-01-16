// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/drug_design/ReactionGrow.cxxtest.hh
/// @brief  test for ReactionGrow Chemistry
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Project Headers
#include <core/types.hh>

#include <protocols/drug_design/ReactionGrow.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/residue_io.hh>

// Utility Headers
#include <utility/file/FileName.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.drug_design.ReactionGrow.cxxtest.hh");

// --------------- Test Class --------------- //

class ReactionGrowTests : public CxxTest::TestSuite {

private:
public:

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	void test_grow() {
		core::chemical::ResidueTypeSetCOP residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD) );
		core::chemical::MutableResidueTypeOP restype( core::chemical::read_topology_file("core/chemical/params/U28.params",residue_set) );

		protocols::drug_design::ReactionGrow grow;
		grow.reaction_file("protocols/drug_design/enamine_rxn.txt" );
		grow.fragment_database("protocols/drug_design/test_fragments1.sdf");

		TS_ASSERT_EQUALS( restype->nheavyatoms(), 26);
		grow.apply(*restype);
		// U28 has hydroxyl, fragment file has a 7 atom aldehyde
		TS_ASSERT_EQUALS( restype->nheavyatoms(), 26+7 );
	}

	//The input molecule should only match the first reactant
	void test_no_first_match() {
		core::chemical::ResidueTypeSetCOP residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD) );
		core::chemical::MutableResidueTypeOP restype( core::chemical::read_topology_file("core/chemical/params/U28.params",residue_set) );

		protocols::drug_design::ReactionGrow grow;
		grow.reaction_file("protocols/drug_design/ester_rxn.txt" );
		grow.fragment_database("protocols/drug_design/test_fragments1.sdf");

		TS_ASSERT_EQUALS( restype->nheavyatoms(), 26);
		grow.apply(*restype);
		// U28 has hydroxyl, fragment file has a 6 atom aldehyde
		TS_ASSERT_EQUALS( restype->nheavyatoms(), 26 );
	}

	// Test if the reaction can't match the ligand at all.
	void test_no_reaction() {
		core::chemical::ResidueTypeSetCOP residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD) );
		core::chemical::MutableResidueTypeOP restype( core::chemical::read_topology_file("core/chemical/params/U13.params",residue_set) );

		protocols::drug_design::ReactionGrow grow;
		grow.reaction_file("protocols/drug_design/ester_rxn2.txt" );
		grow.fragment_database("protocols/drug_design/test_fragments1.sdf");

		TS_ASSERT_EQUALS( restype->nheavyatoms(), 29);
		grow.apply(*restype);
		// U13 has no esterable groups
		TS_ASSERT_EQUALS( restype->nheavyatoms(), 29);
	}

	// Test if the reaction can handle a case where there isn't a matching fragment in the fragment set
	void test_no_fragment() {
		core::chemical::ResidueTypeSetCOP residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD) );
		core::chemical::MutableResidueTypeOP restype( core::chemical::read_topology_file("core/chemical/params/U28.params",residue_set) );

		protocols::drug_design::ReactionGrow grow;
		grow.reaction_file("protocols/drug_design/ester_rxn2.txt" );
		grow.fragment_database("protocols/drug_design/test_fragments2.sdf");

		TS_ASSERT_EQUALS( restype->nheavyatoms(), 26);
		grow.apply(*restype);
		TS_ASSERT_EQUALS( restype->nheavyatoms(), 26 );
	}

	void test_multiple_reaction() {
		// Should find the appropriate reaction, even if there are multiple ones.
		core::chemical::ResidueTypeSetCOP residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD) );
		core::chemical::MutableResidueTypeOP restype( core::chemical::read_topology_file("core/chemical/params/U28.params",residue_set) );

		protocols::drug_design::ReactionGrow grow;
		grow.reaction_file("protocols/drug_design/amide_rxn.txt");
		grow.reaction_file("protocols/drug_design/ester_rxn2.txt", true); // append
		grow.reaction_file("protocols/drug_design/ester_rxn.txt" , true);
		grow.fragment_database("protocols/drug_design/test_fragments1.sdf");

		TS_ASSERT_EQUALS( restype->nheavyatoms(), 26);
		grow.apply(*restype);
		// U28 has hydroxyl, fragment file has a 7 atom aldehyde
		TS_ASSERT_EQUALS( restype->nheavyatoms(), 26+7 );
	}

};
