// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/cyclic_peptide/CycpepSymmetryFilterTests.cxxtest.hh
/// @brief  Unit tests for the CycpepSymmetryFilter.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/cyclic_peptide/CycpepSymmetryFilter.hh>

// Protocols Headers
#include <protocols/cyclic_peptide/DeclareBond.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Protocol Headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("CycpepSymmetryFilterTests");


class CycpepSymmetryFilterTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();
		asymm_peptide_ = core::import_pose::pose_from_file( "protocols/cyclic_peptide/asymm.pdb", false, core::import_pose::PDB_file);
		c2_symm_peptide_ = core::import_pose::pose_from_file( "protocols/cyclic_peptide/c2_symm.pdb", false, core::import_pose::PDB_file);
		c2m_symm_peptide_ = core::import_pose::pose_from_file( "protocols/cyclic_peptide/c2m_symm.pdb", false, core::import_pose::PDB_file);

		protocols::cyclic_peptide::DeclareBond bond;
		bond.set( 8, "C", 1, "N", false );
		bond.apply( *asymm_peptide_ );
		bond.apply( *c2_symm_peptide_ );
		bond.apply( *c2m_symm_peptide_ );

	}

	void tearDown(){

	}

	/// @brief Confirms that the filter picks out the c2-symmetric peptide.
	void test_c2_symm(){
		protocols::cyclic_peptide::CycpepSymmetryFilter filter;
		filter.set_mirror_symm(false);
		filter.set_symm_repeats(2);

		bool const asymm( filter.apply(*asymm_peptide_) );
		bool const c2( filter.apply(*c2_symm_peptide_) );
		bool const c2m( filter.apply(*c2m_symm_peptide_) );

		TR << "c2_symmetry_test:\n\tAsymm\tc2\tc2m";
		TR << "\nExpect:\tFAIL\tPASS\tFAIL";
		TR << "\nActual:\t" << (asymm ? "PASS" : "FAIL") << "\t" << (c2 ? "PASS" : "FAIL") << "\t" << (c2m ? "PASS" : "FAIL") << std::endl;

		TS_ASSERT( !asymm );
		TS_ASSERT( c2 );
		TS_ASSERT( !c2m );
	}

	/// @brief Confirms that the filter picks out the c2m-symmetric peptide.
	void test_c2m_symm(){
		protocols::cyclic_peptide::CycpepSymmetryFilter filter;
		filter.set_mirror_symm(true);
		filter.set_symm_repeats(2);

		bool const asymm( filter.apply(*asymm_peptide_) );
		bool const c2( filter.apply(*c2_symm_peptide_) );
		bool const c2m( filter.apply(*c2m_symm_peptide_) );

		TR << "c2m_symmetry_test:\n\tAsymm\tc2\tc2m";
		TR << "\nExpect:\tFAIL\tFAIL\tPASS";
		TR << "\nActual:\t" << (asymm ? "PASS" : "FAIL") << "\t" << (c2 ? "PASS" : "FAIL") << "\t" << (c2m ? "PASS" : "FAIL") << std::endl;

		TS_ASSERT( !asymm );
		TS_ASSERT( !c2 );
		TS_ASSERT( c2m );
	}

private:

	core::pose::PoseOP asymm_peptide_;
	core::pose::PoseOP c2_symm_peptide_;
	core::pose::PoseOP c2m_symm_peptide_;

};
