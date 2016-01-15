// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/residue_io.cxxtest.hh
/// @brief unit tests for the residue_io file
/// @author Rocco Moretti (rmorettiase@gmail.com)


// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/residue_io.hh>

// Project Headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/MMAtomType.hh>

#include <core/types.hh>

// Platform Headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <string>

#include <test/UTracer.hh>

static basic::Tracer TR("core.chemical.restype_io.cxxtest");

class residue_io_Tests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	void test_graph_out() {
		using namespace core::chemical;

		ChemicalManager * cm(ChemicalManager::get_instance());
		std::string const tag(FA_STANDARD);
		ResidueTypeSetCOP rsd_types = cm->residue_type_set(tag);

		ResidueTypeCOP rsd_ref( rsd_types->name_map("TYR").get_self_ptr() );

		test::UTracer UT("core/chemical/TYR.dot.u");

		write_graphviz( *rsd_ref, UT );
	}

};
