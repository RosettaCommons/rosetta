// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/modifications/Reprotonate.cxxtest.hh
/// @brief  test suite for core::chemical::modifications::Reprotonate.cc
/// @author Rocco Moretti (rmorettiase@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/chemical/modifications/Reprotonate.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/gasteiger/GasteigerAtomTyper.hh>
#include <core/chemical/Element.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/residue_io.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

// Project headers
#include <test/core/init_util.hh>

using namespace core;
using namespace core::chemical;
using namespace core::chemical::gasteiger;

static basic::Tracer TR("core.chemical.modifications.Reprotonate.cxxtest");

// --------------- Test Class --------------- //

class ReprotonateTests : public CxxTest::TestSuite {

public:

	// --------------- Suite-level Fixture --------------- //

	ReprotonateTests() {
	}

	virtual ~ReprotonateTests() {}


	// --------------- Fixtures --------------- //

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// --------------- Test Cases --------------- //

	/// @brief Test that deprotonation works
	void test_deprotonate() {
		using namespace core::chemical;
		ChemicalManager * cm(ChemicalManager::get_instance());
		std::string const tag(FA_STANDARD);
		AtomTypeSetCOP atom_types = cm->atom_type_set(tag);
		ElementSetCOP element_types = cm->element_set("default");
		MMAtomTypeSetCOP mm_atom_types = cm->mm_atom_type_set(tag);
		orbitals::OrbitalTypeSetCOP orbital_types = cm->orbital_type_set(tag);

		core::chemical::MutableResidueTypeOP restype = read_topology_file("core/chemical/modifications/DEP.params",
			atom_types, element_types, mm_atom_types, orbital_types );

		core::chemical::modifications::Reprotonate reprot;

		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("O6")), 1 );  // Carboxylate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("O7")), 0 );  // Carboxylate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("S1")), 1 );  // Thiocarboxylate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("O8")), 0 );  // Thiocarboxylate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("O9")), 1 ); // phosphate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("S2")), 1 ); // phosphate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("O10")), 0 ); // phosphate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("O13")), 1 ); // sulfate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("O12")), 0 ); // sulfate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("S4")), 1 );  // thiophosphate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("O11")), 1 );  // thiophosphate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("S3")), 0 );  // thiophosphate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("O14")), 1 );  // nitrate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("O15")), 0 );  // nitrate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("O17")), 1 );  // enol + vinyl hydroxyl
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("O16")), 1 );  // phenol
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("S6")), 1 );  // thiophenol

		reprot.apply(*restype);

		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("O6")), 0 );  // Carboxylate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("O7")), 0 );  // Carboxylate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("S1")), 0 );  // Thiocarboxylate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("O8")), 0 );  // Thiocarboxylate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("O9")), 0 ); // phosphate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("S2")), 0 ); // phosphate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("O10")), 0 ); // phosphate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("O13")), 0 ); // sulfate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("O12")), 0 ); // sulfate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("S4")), 0 );  // thiophosphate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("O11")), 0 );  // thiophosphate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("S3")), 0 );  // thiophosphate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("O14")), 0 );  // nitrate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("O15")), 0 );  // nitrate
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("O17")), 1 );  // enol + vinyl hydroxyl
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("O16")), 1 );  // phenol
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("S6")), 1 );  // thiophenol

		TS_ASSERT_EQUALS( restype->atom("O6").formal_charge(), -1 );  // Carboxylate
		TS_ASSERT_EQUALS( restype->atom("O7").formal_charge(), 0 );  // Carboxylate
		TS_ASSERT_EQUALS( restype->atom("S1").formal_charge(), -1 );  // Thiocarboxylate
		TS_ASSERT_EQUALS( restype->atom("O8").formal_charge(), 0 );  // Thiocarboxylate
		TS_ASSERT_EQUALS( restype->atom("O9").formal_charge(), -1 ); // phosphate
		TS_ASSERT_EQUALS( restype->atom("S2").formal_charge(), -1 ); // phosphate
		TS_ASSERT_EQUALS( restype->atom("O10").formal_charge(), 0 ); // phosphate
		TS_ASSERT_EQUALS( restype->atom("O13").formal_charge(), -1 ); // sulfate
		TS_ASSERT_EQUALS( restype->atom("O12").formal_charge(), 0 ); // sulfate
		TS_ASSERT_EQUALS( restype->atom("S4").formal_charge(), -1 );  // thiophosphate
		TS_ASSERT_EQUALS( restype->atom("O11").formal_charge(), -1 );  // thiophosphate
		TS_ASSERT_EQUALS( restype->atom("S3").formal_charge(), 0 );  // thiophosphate
		TS_ASSERT_EQUALS( restype->atom("O14").formal_charge(), -1 );  // nitrate
		TS_ASSERT_EQUALS( restype->atom("O15").formal_charge(), 0 );  // nitrate
		TS_ASSERT_EQUALS( restype->atom("O17").formal_charge(), 0 );  // enol + vinyl hydroxyl
		TS_ASSERT_EQUALS( restype->atom("O16").formal_charge(), 0 );  // phenol
		TS_ASSERT_EQUALS( restype->atom("S6").formal_charge(), 0 );  // thiophenol
	}

	/// @brief Test that protonation works
	void test_protonate() {
		using namespace core::chemical;
		ChemicalManager * cm(ChemicalManager::get_instance());
		std::string const tag(FA_STANDARD);
		AtomTypeSetCOP atom_types = cm->atom_type_set(tag);
		ElementSetCOP element_types = cm->element_set("default");
		MMAtomTypeSetCOP mm_atom_types = cm->mm_atom_type_set(tag);
		orbitals::OrbitalTypeSetCOP orbital_types = cm->orbital_type_set(tag);

		core::chemical::MutableResidueTypeOP restype = read_topology_file("core/chemical/modifications/PRO.params",
			atom_types, element_types, mm_atom_types, orbital_types );

		core::chemical::modifications::Reprotonate reprot;

		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("N1")), 2 );  // Primary amine
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("N2")), 1 );  // secondary amine
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("N3")), 0 );  // tertiary amine
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("N4")), 3 );  // Primary amine, already protonated
		TS_ASSERT_EQUALS( restype->atom("N4").formal_charge(), 1 );  // Primary amine, already protonated
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("N5")), 0 );  // quaternary amine
		TS_ASSERT_EQUALS( restype->atom("N5").formal_charge(), 1 );  // quaternary amine
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("N6")), 2 );  // terminal amide
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("N7")), 0 );  // filled amide
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("N8")), 0 );  // cyano
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("N9")), 1 );  // imine
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("N11")), 2 );  // aniline

		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("N10")), 0 );  // aromatic

		reprot.apply(*restype);

		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("N1")), 3 );  // Primary amine
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("N2")), 2 );  // secondary amine
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("N3")), 1 );  // tertiary amine
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("N4")), 3 );  // Primary amine, already protonated
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("N5")), 0 );  // quaternary amine
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("N6")), 2 );  // terminal amide
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("N7")), 0 );  // filled amide
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("N8")), 0 );  // cyano
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("N9")), 1 );  // imine - no change. pKa for conjugate acid is 4-5
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("N11")), 2 );  // aniline - no change. pKa for conjugate acid is 4-5


		TS_ASSERT_EQUALS( restype->atom("N1").formal_charge(), 1 );  // Primary amine
		TS_ASSERT_EQUALS( restype->atom("N2").formal_charge(), 1 );  // secondary amine
		TS_ASSERT_EQUALS( restype->atom("N3").formal_charge(), 1 );  // tertiary amine
		TS_ASSERT_EQUALS( restype->atom("N4").formal_charge(), 1 );  // Primary amine, already protonated
		TS_ASSERT_EQUALS( restype->atom("N5").formal_charge(), 1 );  // quaternary amine
		TS_ASSERT_EQUALS( restype->atom("N6").formal_charge(), 0 );  // terminal amide
		TS_ASSERT_EQUALS( restype->atom("N7").formal_charge(), 0 );  // filled amide
		TS_ASSERT_EQUALS( restype->atom("N8").formal_charge(), 0 );  // cyano
		TS_ASSERT_EQUALS( restype->atom("N9").formal_charge(), 0 );  // imine
		TS_ASSERT_EQUALS( restype->atom("N11").formal_charge(), 0 );  // aniline

		// Depending on system, protonation may be possible. However, we'll punt and leave aromatics as-is
		// (Anyway, pyridine conjugate acid is pKa is ~5.3)
		TS_ASSERT_EQUALS( restype->number_bonded_hydrogens(restype->atom_vertex("N10")), 0 );  // aromatic
		TS_ASSERT_EQUALS( restype->atom("N10").formal_charge(), 0 );  // aromatic
	}
};


