// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/ResidueTypeTests.cxxtest.hh
/// @brief unit tests for ResidueType
/// @author Matthew O'Meara


// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/chemical/ResidueType.hh>

// Project Headers
// AUTO-REMOVED #include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/chemical/ElementSet.hh>
// AUTO-REMOVED #include <core/chemical/MMAtomTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/orbitals/OrbitalTypeSet.hh>

// Platform Headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <string>
#include <ostream>

//Auto Headers


using std::endl;
using std::string;
using basic::Tracer;
using core::chemical::AtomTypeSetCAP;
using core::chemical::ChemicalManager;
using core::chemical::ElementSetCAP;
using core::chemical::MMAtomTypeSetCAP;
using core::chemical::orbitals::OrbitalTypeSetCAP;
using core::chemical::ResidueType;
using core::chemical::FA_STANDARD;
using utility::vector1;

static Tracer TR("core.chemical.ResidueTypeTests.cxxtest");

class ResidueTypeTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	void test_residue_type_initialization() {
		ChemicalManager * cm(ChemicalManager::get_instance());

		string const tag(FA_STANDARD);

		// minirosett_database/chemical/atom_type_sets/<tag>
		AtomTypeSetCAP atom_types = cm->atom_type_set(tag);

		// minirosetta_database/chemical/element_sets/<tag>
		ElementSetCAP element_types = cm->element_set(tag);

		// minirosetta_database/chemical/mm_atom_type_sets/<tag>
		MMAtomTypeSetCAP mm_atom_types = cm->mm_atom_type_set(tag);

		// minirosetta_database/chemical/orbital_type_sets/<tag>
		OrbitalTypeSetCAP orbital_types = cm->orbital_type_set(tag);

		ResidueType rsd( atom_types, element_types, mm_atom_types, orbital_types);

		string p("property1");
		rsd.add_property(p);
		vector1<string > const & properties(rsd.properties());
 		TS_ASSERT(!properties[1].compare(p));

	}

};
