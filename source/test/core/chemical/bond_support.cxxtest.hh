// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/chemical/bond_support.cxxtest.hh
/// @brief unit tests for the bond_support file
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/bond_support.hh>
#include <core/chemical/Atom.hh>
#include <core/chemical/Bond.hh>

// Project Headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/residue_io.hh>
#include <core/types.hh>

// Platform Headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <string>

static basic::Tracer TR("core.chemical.bond_support.cxxtest");

class bond_support_Tests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	void test_find_bonds_in_rings() {
		using namespace core::chemical;
		ChemicalManager * cm(ChemicalManager::get_instance());
		std::string const tag(FA_STANDARD);
		AtomTypeSetCOP atom_types = cm->atom_type_set(tag);
		ElementSetCOP element_types = cm->element_set("default");
		MMAtomTypeSetCOP mm_atom_types = cm->mm_atom_type_set(tag);

		ResidueType res( atom_types, element_types, mm_atom_types, NULL );

		TR << "Testing Biphenyl" << std::endl;
		//Biphenyl - we want the ring bonds but not the connecting or hygrogen bond to be rings
		//We shouldn't need much more than the atom layout.
		res.add_atom("CA1");
		res.add_atom("CA2");
		res.add_atom("CA3");
		res.add_atom("CA4");
		res.add_atom("CA5");
		res.add_atom("CA6");
		res.add_atom("CB1");
		res.add_atom("CB2");
		res.add_atom("CB3");
		res.add_atom("CB4");
		res.add_atom("CB5");
		res.add_atom("CB6");
		res.add_atom("HA2");
		res.add_atom("HA3");
		res.add_atom("HA4");
		res.add_atom("HA5");
		res.add_atom("HA6");
		res.add_atom("HB2");
		res.add_atom("HB3");
		res.add_atom("HB4");
		res.add_atom("HB5");
		res.add_atom("HB6");
		// Ring 1
		res.add_bond("CA1","CA2");
		res.add_bond("CA2","CA3", DoubleBond);
		res.add_bond("CA3","CA4");
		res.add_bond("CA4","CA5", DoubleBond);
		res.add_bond("CA5","CA6");
		res.add_bond("CA6","CA1", DoubleBond);
		//Bridge
		res.add_bond("CA1","CB1");
		//Ring 2
		res.add_bond("CB1","CB2");
		res.add_bond("CB2","CB3", DoubleBond);
		res.add_bond("CB3","CB4");
		res.add_bond("CB4","CB5", DoubleBond);
		res.add_bond("CB5","CB6");
		res.add_bond("CB6","CB1", DoubleBond);
		//Hydrogens
		res.add_bond("CA2","HA2");
		res.add_bond("CA3","HA3");
		res.add_bond("CA4","HA4");
		res.add_bond("CA5","HA5");
		res.add_bond("CA6","HA6");
		res.add_bond("CB2","HB2");
		res.add_bond("CB3","HB3");
		res.add_bond("CB4","HB4");
		res.add_bond("CB5","HB5");
		res.add_bond("CB6","HB6");

		find_bonds_in_rings(res);

		//Ring 1
		TS_ASSERT( res.bond("CA1","CA2").ringness() == BondInRing );
		TS_ASSERT( res.bond("CA2","CA3").ringness() == BondInRing );
		TS_ASSERT( res.bond("CA3","CA4").ringness() == BondInRing );
		TS_ASSERT( res.bond("CA4","CA5").ringness() == BondInRing );
		TS_ASSERT( res.bond("CA5","CA6").ringness() == BondInRing );
		TS_ASSERT( res.bond("CA6","CA1").ringness() == BondInRing );
		//Bridge
		TS_ASSERT( res.bond("CA1","CB1").ringness() == BondNotInRing );
		//Ring 2
		TS_ASSERT( res.bond("CB1","CB2").ringness() == BondInRing );
		TS_ASSERT( res.bond("CB2","CB3").ringness() == BondInRing );
		TS_ASSERT( res.bond("CB3","CB4").ringness() == BondInRing );
		TS_ASSERT( res.bond("CB4","CB5").ringness() == BondInRing );
		TS_ASSERT( res.bond("CB5","CB6").ringness() == BondInRing );
		TS_ASSERT( res.bond("CB6","CB1").ringness() == BondInRing );
		//Hydrogens
		TS_ASSERT( res.bond("CA2","HA2").ringness() == BondNotInRing );
		TS_ASSERT( res.bond("CA3","HA3").ringness() == BondNotInRing );
		TS_ASSERT( res.bond("CA4","HA4").ringness() == BondNotInRing );
		TS_ASSERT( res.bond("CA5","HA5").ringness() == BondNotInRing );
		TS_ASSERT( res.bond("CA6","HA6").ringness() == BondNotInRing );
		TS_ASSERT( res.bond("CB2","HB2").ringness() == BondNotInRing );
		TS_ASSERT( res.bond("CB3","HB3").ringness() == BondNotInRing );
		TS_ASSERT( res.bond("CB4","HB4").ringness() == BondNotInRing );
		TS_ASSERT( res.bond("CB5","HB5").ringness() == BondNotInRing );
		TS_ASSERT( res.bond("CB6","HB6").ringness() == BondNotInRing );
	}

	void test_find_bonds_in_rings2() {
		using namespace core::chemical;
		ChemicalManager * cm(ChemicalManager::get_instance());
		ResidueTypeSetCOP restypeset = cm->residue_type_set(FA_STANDARD);

		ResidueTypeOP res( read_topology_file( "core/chemical/params/U26.params", restypeset ) );

		TR << "Testing U26" << std::endl;

		find_bonds_in_rings(*res);

		//Bridging bond.
		TS_ASSERT( res->bond("C13","C14").ringness() == BondNotInRing );
		// Exocyclic guanidinium
		TS_ASSERT( res->bond("C15","N4").ringness() == BondNotInRing );
		TS_ASSERT( res->bond("C16","N4").ringness() == BondNotInRing );
		TS_ASSERT( res->bond("C16","N2").ringness() == BondNotInRing );
		TS_ASSERT( res->bond("C16","N3").ringness() == BondNotInRing );
		// Internal nitrogen
		TS_ASSERT( res->bond("C15","N1").ringness() == BondInRing );
		TS_ASSERT( res->bond("C10","N1").ringness() == BondInRing );
		// Fused ring bond
		TS_ASSERT( res->bond("C11","C12").ringness() == BondInRing );
	}

};
