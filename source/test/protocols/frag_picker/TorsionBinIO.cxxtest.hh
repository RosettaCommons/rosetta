// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/TorsionBinIO.cxxtest.hh
/// @brief
/// @author James Thompson

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>
#include <protocols/frag_picker/TorsionBinIO.hh>
#include <utility/io/izstream.hh>

//Auto Headers
#include <utility/vector1.hh>


static basic::Tracer tr("protocols.frag_picker.TorsionBinIO.cxxtest");
//MY_TRACERS("core.fragment.ConstantLengthFragments.cxxtest")

class TorsionBinIOTests : public CxxTest::TestSuite {
public:
	TorsionBinIOTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// tests
	void test_read_file() {
		protocols::frag_picker::TorsionBinIO tb_io;
		tb_io.read( utility::io::izstream( "protocols/frag_picker/t286_.6.csts" ) );
		core::Real const DELTA( 1e-5 );

		TS_ASSERT( tb_io.nrows() == 200 );
		TS_ASSERT_DELTA( tb_io.prof_row(1)  [5], 0.006,    DELTA );
		TS_ASSERT_DELTA( tb_io.prof_row(20) [3], 0.057022, DELTA );
		TS_ASSERT_DELTA( tb_io.prof_row(200)[2], 0.782440, DELTA );
	}

	// Shared finalization goes here.
	void tearDown() {}

private:

}; // TorsionBinIOTests
