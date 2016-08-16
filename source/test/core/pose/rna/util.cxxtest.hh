// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/pose/rna/util.cxxtest.hh
/// @brief  unit tests for core::pose::rna::util.hh/cc functions
/// @author Andy Watkins

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <core/types.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/util.hh>


#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>

#include <test/util/pose_funcs.hh>

//Auto Headers
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>


class PoseRNAUtilTests : public CxxTest::TestSuite {


public: // setup


	typedef core::Size Size;
	typedef core::pose::PDBInfo PDBInfo;
	typedef core::pose::PDBInfoOP PDBInfoOP;
	typedef core::pose::Pose Pose;
	typedef core::pose::PoseOP PoseOP;


public: //setup


	PoseRNAUtilTests() {};


	// Shared initialization.
	void setUp() {
		core_init();
	}


	// Shared finalization.
	void tearDown() {}

public: // re-used methods

	void test_remove_nonnaturals() {
		std::string seq = "acguacguX[2MA]X[1MA]acguacgu";
		std::string cleanseq;
		std::map< Size, std::string > special_res;
		core::pose::rna::remove_and_store_bracketed( seq, cleanseq, special_res);

		TS_ASSERT( cleanseq == "acguacguXXacguacgu" );
		TS_ASSERT( special_res[ 8 ] == "2MA" );
		TS_ASSERT( special_res[ 9 ] == "1MA" );

		// cleaning a clean sequence should be a noop
		std::string seq2 = "acguacgu";
		TS_ASSERT( seq2 == core::pose::rna::remove_bracketed( seq2 ) );

	}

}; // class PoseRNAUtilTests
