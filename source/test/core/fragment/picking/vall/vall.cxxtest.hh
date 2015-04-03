// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/fragment/picking_old/vall/vall.cxxtest.hh
/// @brief  unit tests for various Vall library functionalities
/// @remarks this file will likely get broken up into multiple unit tests as
///  development proceeds and the picking_old framework is moved into core
/// @author

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <core/types.hh>
#include <core/fragment/picking_old/vall/vall_io.hh>
#include <core/fragment/picking_old/vall/VallSection.hh>
#include <core/fragment/picking_old/vall/VallLibrary.hh>
#include <core/fragment/picking_old/vall/VallLibrarian.hh>
#include <core/fragment/picking_old/vall/eval/IdentityEval.hh>
#include <core/fragment/picking_old/vall/gen/LengthGen.hh>
#include <core/fragment/picking_old/vall/gen/SecStructGen.hh>

#include <string>

//Auto Headers
#include <platform/types.hh>
#include <core/chemical/AA.hh>
#include <core/fragment/BBTorsionSRFD.fwd.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/FragData.fwd.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/fragment/SecstructSRFD.hh>
#include <core/fragment/SingleResidueFragData.fwd.hh>
#include <core/fragment/SingleResidueFragData.hh>
#include <core/fragment/picking_old/concepts/Book.fwd.hh>
#include <core/fragment/picking_old/concepts/Book.hh>
#include <core/fragment/picking_old/concepts/Extent.fwd.hh>
#include <core/fragment/picking_old/concepts/Extent.hh>
#include <core/fragment/picking_old/concepts/Librarian.fwd.hh>
#include <core/fragment/picking_old/concepts/Librarian.hh>
#include <core/fragment/picking_old/concepts/Library.fwd.hh>
#include <core/fragment/picking_old/concepts/Library.hh>
#include <core/fragment/picking_old/vall/VallLibrarian.fwd.hh>
#include <core/fragment/picking_old/vall/VallLibrary.fwd.hh>
#include <core/fragment/picking_old/vall/VallResidue.fwd.hh>
#include <core/fragment/picking_old/vall/VallResidue.hh>
#include <core/fragment/picking_old/vall/VallSection.fwd.hh>
#include <core/fragment/picking_old/vall/eval/IdentityEval.fwd.hh>
#include <core/fragment/picking_old/vall/eval/VallFragmentEval.fwd.hh>
#include <core/fragment/picking_old/vall/eval/VallFragmentEval.hh>
#include <core/fragment/picking_old/vall/gen/LengthGen.fwd.hh>
#include <core/fragment/picking_old/vall/gen/SecStructGen.fwd.hh>
#include <core/fragment/picking_old/vall/gen/VallFragmentGen.fwd.hh>
#include <core/fragment/picking_old/vall/gen/VallFragmentGen.hh>
#include <core/fragment/picking_old/vall/scores/VallFragmentScore.fwd.hh>
#include <core/fragment/picking_old/vall/scores/VallFragmentScore.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/fix_boinc_read.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/xyzVector.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <iosfwd>
#include <iostream>
#include <iterator>
#include <limits>
#include <ostream>
#include <sstream>
#include <vector>
#include <basic/Tracer.fwd.hh>
#include <basic/Tracer.hh>


class VallTests : public CxxTest::TestSuite
{


public: //setup


	VallTests() {};


	// Shared initialization.
	void setUp() {
		core_init();
	}


	// Shared finalization.
	void tearDown() {
	}


public: // tests
	/// @brief test standard read function
	void test_vall_library_from_file() {
		using namespace core::fragment::picking_old::vall;

		// read the test library
		VallLibrary library;
		vall_library_from_file( "core/fragment/picking/vall/vall_test.dat", library );

		// test library sizes
		TS_ASSERT_EQUALS( library.size(), 3 );
		TS_ASSERT_EQUALS( library.n_residues(), 16 );

		// test per-book sizes
		VallLibrary::BookConstIterator iter = library.begin();
		TS_ASSERT_EQUALS( iter->size(), 4 );
		++iter;
		TS_ASSERT_EQUALS( iter->size(), 5 );
		++iter;
		TS_ASSERT_EQUALS( iter->size(), 7 );
		++iter;
		TS_ASSERT( iter == library.end() );

		// test content of a page
		iter = ++library.begin(); // 2nd book
		VallSection::PageConstIterator p = ++( iter->begin() ); // 2nd page
		TS_ASSERT_EQUALS( p->id(), "8tlnE" );
		TS_ASSERT_EQUALS( p->aa(), 'A' );
		TS_ASSERT_EQUALS( p->ss(), 'H' );
		TS_ASSERT_EQUALS( p->resi(), 912 );
		TS_ASSERT_EQUALS( p->x(), 48.61 );
		TS_ASSERT_EQUALS( p->y(), 17.34 );
		TS_ASSERT_EQUALS( p->z(), -5.58 );
		TS_ASSERT_EQUALS( p->phi(), -66.477 );
		TS_ASSERT_EQUALS( p->psi(), -25.529 );
		TS_ASSERT_EQUALS( p->omega(), 178.457 );
		TS_ASSERT_EQUALS( p->profile()[ p->profile().size() ], 0.008 ); // last profile value

	}


	/// @brief test cataloging function
	void test_catalog() {
		using namespace core::fragment::picking_old::vall;
		using namespace core::fragment::picking_old::vall::gen;
		using namespace core::fragment::picking_old::vall::eval;

		// read the test library
		VallLibrary library;
		vall_library_from_file( "core/fragment/picking/vall/vall_test.dat", library );

		VallLibrarian librarian;
		librarian.add_fragment_gen( VallFragmentGenCOP( new gen::LengthGen( 3 ) ) ); // 3-mers
		librarian.add_fragment_eval( VallFragmentEvalCOP( new eval::IdentityEval( "HLE", "...", 1.0, 0.0, true ) ) );

		// catalog fragment library sorting via bookmark scores
		librarian.catalog( library );

		TS_ASSERT_EQUALS( librarian.n_scores(), 10 );
	}


	/// @brief test SecStructGen
	void test_SecStructGen() {
		using namespace core::fragment::picking_old::vall;
		using namespace core::fragment::picking_old::vall::gen;
		using namespace core::fragment::picking_old::vall::eval;

		// read the test library
		VallLibrary library;
		vall_library_from_file( "core/fragment/picking/vall/vall_test.dat", library );

		VallLibrarian librarian;
		librarian.add_fragment_gen( VallFragmentGenCOP( new gen::SecStructGen( "EEH" ) ) ); // 3-mer
		librarian.add_fragment_eval( VallFragmentEvalCOP( new eval::IdentityEval( "EEH", "...", 1.0, 0.0, true ) ) );

		// catalog fragment library sorting via bookmark scores
		librarian.catalog( library );

		// there should be only 1 score since only a single EEH 3-mer exists in
		// the vall_test.dat file
		TS_ASSERT_EQUALS( librarian.n_scores(), 1 );
	}


};
