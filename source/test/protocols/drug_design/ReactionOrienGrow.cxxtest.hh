// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/drug_design/ReactionOrienGrow.cxxtest.hh
/// @brief  test for ReactionOrienGrow Chemistry
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Project Headers
#include <core/types.hh>

#include <protocols/drug_design/ReactionOrienGrow.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/residue_io.hh>

// Utility Headers
#include <utility/file/FileName.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.drug_design.ReactionOrienGrow.cxxtest.hh");

// --------------- Test Class --------------- //

class ReactionOrienGrowTests : public CxxTest::TestSuite {

private:
public:

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	void test_grow() {
		core::chemical::ResidueTypeSetCOP residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD) );
		core::chemical::ResidueTypeOP restype( core::chemical::read_topology_file("core/chemical/params/U28.params",residue_set) );

		protocols::drug_design::ReactionOrienGrow grow;
		grow.load_reactions("protocols/drug_design", "enamine_rxn.txt" );

		//TS_ASSERT_EQUALS( restype->nheavyatoms(), 26);
		grow.apply(*restype);
		// Reagent 1 have 11 heavy atoms + reagent 2 have 15 heavy atoms. 
		// Two atoms lost during the reaction.
		TS_ASSERT_EQUALS( restype->nheavyatoms(), 11+15-2 );
	}

};

