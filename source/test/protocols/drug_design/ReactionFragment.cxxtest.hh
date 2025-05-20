// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/drug_design/ReactionFragment.cxxtest.hh
/// @brief  test for ReactionFragment Chemistry
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Project Headers
#include <core/types.hh>

#include <protocols/drug_design/ReactionFragment.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/residue_io.hh>

// Utility Headers
#include <utility/file/FileName.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.drug_design.ReactionFragment.cxxtest.hh");

// --------------- Test Class --------------- //

class ReactionFragmentTests : public CxxTest::TestSuite {

private:
public:

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	void test_default() {
		core::chemical::ResidueTypeSetCOP residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD) );
		{
			core::chemical::MutableResidueTypeOP restype( core::chemical::read_topology_file("core/chemical/params/U27.params",residue_set) );
			core::chemical::NameVDMapping orig_name_map( *restype );

			protocols::drug_design::ReactionFragment fragment;
			fragment.reaction_file("protocols/drug_design/ester_rxn.txt" );

			TS_ASSERT_EQUALS( restype->nheavyatoms(), 18);
			fragment.apply(*restype);
			core::chemical::NameVDMapping mapping( orig_name_map.downstream_combine( fragment.get_mapping()) );

			TS_ASSERT_EQUALS( restype->nheavyatoms(), 12);
			//TS_ASSERT_DIFFERS( mapping[ " O1 " ], core::chemical::MutableResidueType::null_vertex ); // Atoms involved are currently ignored
			//TS_ASSERT_EQUALS( mapping[ " O2 " ], core::chemical::MutableResidueType::null_vertex ); // Atoms involved are currently ignored
			TS_ASSERT_DIFFERS( mapping[ " N1 " ], core::chemical::MutableResidueType::null_vertex );
			TS_ASSERT_EQUALS( mapping[ " C13" ], core::chemical::MutableResidueType::null_vertex );
		}
		{
			core::chemical::MutableResidueTypeOP restype( core::chemical::read_topology_file("core/chemical/params/U27.params",residue_set) );
			core::chemical::NameVDMapping orig_name_map( *restype );

			protocols::drug_design::ReactionFragment fragment;
			fragment.reaction_file("protocols/drug_design/ester_rxn2.txt" );

			TS_ASSERT_EQUALS( restype->nheavyatoms(), 18);
			fragment.apply(*restype);
			core::chemical::NameVDMapping mapping( orig_name_map.downstream_combine( fragment.get_mapping()) );

			TS_ASSERT_EQUALS( restype->nheavyatoms(), 6);
			//TS_ASSERT_EQUALS( mapping[ " O1 " ], core::chemical::MutableResidueType::null_vertex ); // Atoms involved are currently ignored
			//TS_ASSERT_DIFFERS( mapping[ " O2 " ], core::chemical::MutableResidueType::null_vertex ); // Atoms involved are currently ignored
			TS_ASSERT_EQUALS( mapping[ " N1 " ], core::chemical::MutableResidueType::null_vertex );
			TS_ASSERT_DIFFERS( mapping[ " C13" ], core::chemical::MutableResidueType::null_vertex );
		}
	}

	void test_larger() {
		core::chemical::ResidueTypeSetCOP residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD) );
		{
			core::chemical::MutableResidueTypeOP restype( core::chemical::read_topology_file("core/chemical/params/U13.params",residue_set) );
			core::chemical::NameVDMapping orig_name_map( *restype );

			protocols::drug_design::ReactionFragment fragment;
			fragment.reaction_file("protocols/drug_design/amide_rxn.txt" );
			fragment.keep_bigger(true);

			TS_ASSERT_EQUALS( restype->nheavyatoms(), 29);
			fragment.apply(*restype);
			core::chemical::NameVDMapping mapping( orig_name_map.downstream_combine( fragment.get_mapping()) );

			TS_ASSERT_EQUALS( restype->nheavyatoms(), 15); // The carbonyl one.
			//TS_ASSERT_DIFFERS( mapping[ " O1 " ], core::chemical::MutableResidueType::null_vertex ); // Atoms involved are currently ignored
			//TS_ASSERT_EQUALS( mapping[ " N2 " ], core::chemical::MutableResidueType::null_vertex ); // Atoms involved are currently ignored
			TS_ASSERT_DIFFERS( mapping[ " N1 " ], core::chemical::MutableResidueType::null_vertex );
			TS_ASSERT_EQUALS( mapping[ " C21" ], core::chemical::MutableResidueType::null_vertex );
		}
		{
			core::chemical::MutableResidueTypeOP restype( core::chemical::read_topology_file("core/chemical/params/U27.params",residue_set) );
			core::chemical::NameVDMapping orig_name_map( *restype );

			protocols::drug_design::ReactionFragment fragment;
			fragment.reaction_file("protocols/drug_design/ester_rxn2.txt" );
			fragment.keep_bigger(true);

			TS_ASSERT_EQUALS( restype->nheavyatoms(), 18);
			fragment.apply(*restype);
			core::chemical::NameVDMapping mapping( orig_name_map.downstream_combine( fragment.get_mapping()) );

			TS_ASSERT_EQUALS( restype->nheavyatoms(), 12);
			//TS_ASSERT_DIFFERS( mapping[ " O1 " ], core::chemical::MutableResidueType::null_vertex ); // Atoms involved are currently ignored
			//TS_ASSERT_EQUALS( mapping[ " O2 " ], core::chemical::MutableResidueType::null_vertex ); // Atoms involved are currently ignored
			TS_ASSERT_DIFFERS( mapping[ " N1 " ], core::chemical::MutableResidueType::null_vertex );
			TS_ASSERT_EQUALS( mapping[ " C13" ], core::chemical::MutableResidueType::null_vertex );
		}
	}

	void test_keep_atom() {
		core::chemical::ResidueTypeSetCOP residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD) );
		core::chemical::MutableResidueTypeOP restype( core::chemical::read_topology_file("core/chemical/params/U13.params",residue_set) );
		core::chemical::NameVDMapping orig_name_map( *restype );

		protocols::drug_design::ReactionFragment fragment;
		fragment.reaction_file("protocols/drug_design/amide_rxn.txt" );
		fragment.keep_bigger(false);
		//fragment.keep_atom("N4"); // Positivly charges - will be lost in conversion
		fragment.keep_atom("C21");

		TS_ASSERT_EQUALS( restype->nheavyatoms(), 29);
		TS_ASSERT( restype->has("C21"));
		fragment.apply(*restype);
		core::chemical::NameVDMapping mapping( orig_name_map.downstream_combine( fragment.get_mapping()) );

		TS_ASSERT_EQUALS( restype->nheavyatoms(), 14);
		TS_ASSERT_DIFFERS( mapping[ " C21" ], core::chemical::MutableResidueType::null_vertex );
		TS_ASSERT_EQUALS( mapping[ " N1 " ], core::chemical::MutableResidueType::null_vertex );
		//TS_ASSERT_DIFFERS( mapping[ " N2 " ], core::chemical::MutableResidueType::null_vertex ); // Atoms involved are currently ignored
		//TS_ASSERT_EQUALS( mapping[ " O1 " ], core::chemical::MutableResidueType::null_vertex ); // Atoms involved are currently ignored
	}

	void test_bad_reaction() {
		// A reaction which doesn't apply to the residue.
		core::chemical::ResidueTypeSetCOP residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD) );
		core::chemical::MutableResidueTypeOP restype( core::chemical::read_topology_file("core/chemical/params/U13.params",residue_set) );
		core::chemical::NameVDMapping orig_name_map( *restype );

		protocols::drug_design::ReactionFragment fragment;
		fragment.reaction_file("protocols/drug_design/ester_rxn.txt");
		fragment.keep_bigger(false);
		//fragment.keep_atom("N4"); // Positivly charges - will be lost in conversion
		fragment.keep_atom("C21");

		TS_ASSERT_EQUALS( restype->nheavyatoms(), 29);
		TS_ASSERT( restype->has("C21"));
		fragment.apply(*restype); // Should be a no-op
		core::chemical::NameVDMapping mapping( orig_name_map.downstream_combine( fragment.get_mapping()) );

		TS_ASSERT_EQUALS( restype->nheavyatoms(), 29);
		TS_ASSERT_DIFFERS( mapping[ " C21" ], core::chemical::MutableResidueType::null_vertex );
		TS_ASSERT_DIFFERS( mapping[ " N2 " ], core::chemical::MutableResidueType::null_vertex );
		TS_ASSERT_DIFFERS( mapping[ " O1 " ], core::chemical::MutableResidueType::null_vertex );
		TS_ASSERT_DIFFERS( mapping[ " N1 " ], core::chemical::MutableResidueType::null_vertex );
	}

	void test_multiple_reaction() {
		// Should find the appropriate reaction, even if there are multiple ones.
		core::chemical::ResidueTypeSetCOP residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD) );
		core::chemical::MutableResidueTypeOP restype( core::chemical::read_topology_file("core/chemical/params/U13.params",residue_set) );
		core::chemical::NameVDMapping orig_name_map( *restype );

		protocols::drug_design::ReactionFragment fragment;
		fragment.reaction_file("protocols/drug_design/ester_rxn.txt", /*append=*/ true );
		fragment.reaction_file("protocols/drug_design/amide_rxn.txt", /*append=*/ true  );
		fragment.keep_bigger(false);
		//fragment.keep_atom("N4"); // Positivly charges - will be lost in conversion
		fragment.keep_atom("C21");

		TS_ASSERT_EQUALS( restype->nheavyatoms(), 29);
		TS_ASSERT( restype->has("C21"));
		fragment.apply(*restype);
		core::chemical::NameVDMapping mapping( orig_name_map.downstream_combine( fragment.get_mapping()) );

		TS_ASSERT_EQUALS( restype->nheavyatoms(), 14);
		TS_ASSERT_DIFFERS( mapping[ " C21" ], core::chemical::MutableResidueType::null_vertex );
		TS_ASSERT_EQUALS( mapping[ " N1 " ], core::chemical::MutableResidueType::null_vertex );
		//TS_ASSERT_DIFFERS( mapping[ " N2 " ], core::chemical::MutableResidueType::null_vertex ); // Atoms involved are currently ignored
		//TS_ASSERT_EQUALS( mapping[ " O1 " ], core::chemical::MutableResidueType::null_vertex ); // Atoms involved are currently ignored
	}

};
