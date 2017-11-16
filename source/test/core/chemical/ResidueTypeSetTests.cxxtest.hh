// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/ResidueTypeSetTests.cxxtest.hh
/// @brief unit test for ResidueTypeSet, not really complete yet
/// @author Florian Richter, jan 11,


// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/GlobalResidueTypeSet.hh>
#include <core/chemical/ResidueTypeFinder.hh>

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

using std::endl;
using std::string;

static basic::Tracer TR("core.chemical.ResidueTypeSetTests.cxxtest");

class ResidueTypeSetTests : public CxxTest::TestSuite {

public:
	core::Real delta_percent;

	void setUp() {
		core_init();
		delta_percent=0.0001;
	}

	void tearDown() {}

	void test_residue_type_sets() {
		using namespace core::chemical;
		using namespace ObjexxFCL::format;

		int width = 15;

		GlobalResidueTypeSetOP rs( new GlobalResidueTypeSet( FA_STANDARD, basic::database::full_name( "chemical/residue_type_sets/"+FA_STANDARD+"/" ) ) );

		TR << A(width,"ResidueTypeSet") << A(width,"NumBaseResTypes") << endl;
		TR << A(width, FA_STANDARD) << I(width,rs->base_residue_types().size()) << endl;
		TR << A(width,"ResidueTypeSet") << A(width,"NumCustomResTypes") << endl;
		TR << A(width, FA_STANDARD) << I(width,rs->unpatchable_residue_types().size()) << endl;

		ResidueType const & serine = rs->name_map( "SER" );
		TS_ASSERT_DELTA(serine.mass(), 87.0900, delta_percent);

	}

	/// @brief Tests the ability of the ResidueTypeSet to give me the mirror-image ResidueType to a given type.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void test_get_mirror_type() {
		using namespace core::chemical;

		GlobalResidueTypeSetOP restypeset( new GlobalResidueTypeSet( FA_STANDARD, basic::database::full_name( "chemical/residue_type_sets/"+FA_STANDARD+"/" ) ) );

		//Test types:
		ResidueTypeCOP lcys_nterm_cterm( restypeset->name_mapOP("CYS:NtermProteinFull:CtermProteinFull") );
		ResidueTypeCOP dile_nacetyl( restypeset->name_mapOP("DILE:AcetylatedNtermProteinFull") );
		ResidueTypeCOP plain_c53( restypeset->name_mapOP("C53") );

		//Return types:
		ResidueTypeCOP mirrored_lcys_nterm_cterm( restypeset->get_mirrored_type( lcys_nterm_cterm ) );
		ResidueTypeCOP mirrored_dile_nacetyl( restypeset->get_mirrored_type( dile_nacetyl ) );
		ResidueTypeCOP mirrored_plain_c53( restypeset->get_mirrored_type( plain_c53 ) );

		TS_ASSERT_EQUALS( mirrored_lcys_nterm_cterm->name3(), "DCS" );
		TS_ASSERT_EQUALS( mirrored_lcys_nterm_cterm->name(), "DCYS:CtermProteinFull:NtermProteinFull" );
		TS_ASSERT_EQUALS( mirrored_lcys_nterm_cterm->base_name(), "DCYS" );

		TS_ASSERT_EQUALS( mirrored_dile_nacetyl->name3(), "ILE" );
		TS_ASSERT_EQUALS( mirrored_dile_nacetyl->name(), "ILE:AcetylatedNtermProteinFull" );
		TS_ASSERT_EQUALS( mirrored_dile_nacetyl->base_name(), "ILE" );

		TS_ASSERT_EQUALS( mirrored_plain_c53->name3(), "C53" ); //Currently, the D-patches don't update the name3 for noncanonicals.  Fix this if/when they do.
		TS_ASSERT_EQUALS( mirrored_plain_c53->name(), "DC53" );
		TS_ASSERT_EQUALS( mirrored_plain_c53->base_name(), "DC53" );

	}

};
