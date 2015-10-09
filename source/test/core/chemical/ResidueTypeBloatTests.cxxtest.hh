// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/ResidueTypeBloat.cxxtest.hh
/// @brief  Measure memory contribution for Residue Types
/// @author Matthew O'Meara


// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/chemical/ResidueTypeSet.hh>

// Project Headers
#include <core/chemical/ChemicalManager.hh>

// Platform Headers
#include <basic/Tracer.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ Headers
#include <string>
#include <ostream>

//Auto Headers
#include <utility/vector1.hh>


using std::endl;
using std::string;

static basic::Tracer TR("core.chemical.ResidueTypeBloatTests.cxxtest");

class ResidueTypeBloatTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	void test_number_of_residue_types() {
		using namespace core::chemical;
		using namespace ObjexxFCL::format;

		int width = 15;
		TR << A(width,"ResidueTypeSet") << A(width,"NumBaseResTypes") << A(width,"NumPatches") << endl;

		string rss;
		ResidueTypeSetCOP rs;

		rss = FA_STANDARD;
		rs = ChemicalManager::get_instance()->residue_type_set(rss );
		TR << A(width, rss)
			<< I(width,rs->base_residue_types().size())
			<< I(width,rs->patches().size())
			<< endl;

		rss = CENTROID;
		rs = ChemicalManager::get_instance()->residue_type_set(rss );
		TR << A(width, rss)
			<< I(width,rs->base_residue_types().size())
			<< I(width,rs->patches().size())
			<< endl;

		// Is this broken?
		//rss = COARSE_TWO_BEAD;
		//rs = ChemicalManager::get_instance()->residue_type_set(rss );
		//TR << A(width, rss) << I(width,rs->residue_types_DO_NOT_USE().size()) << endl;

		// Is this broken?
		//rss = HYBRID_FA_STANDARD_CENTROID;
		//rs = ChemicalManager::get_instance()->residue_type_set(rss );
		//TR << A(width, rss) << I(width,rs->residue_types_DO_NOT_USE().size()) << endl;


	}


};
