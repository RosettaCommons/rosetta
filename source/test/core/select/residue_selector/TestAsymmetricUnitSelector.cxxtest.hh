// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/select/residue_selector/TestMasterSubunitSelector.cxxtest.hh
/// @brief  Unit test for MasterSubunitSelector
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/select/residue_selector/AsymmetricUnitSelector.hh>

// Project Headers


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <test/util/pose_funcs.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("TestAsymmetricUnitSelector");


class TestAsymmetricUnitSelector : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();

	}

	void test_select_master_subunit_selector(){
		using namespace core::conformation::symmetry;
		using namespace core::select::residue_selector;

		core::pose::Pose start_pose;
		core::import_pose::pose_from_file( start_pose, "core/scoring/symmetry/test_in.pdb" , core::import_pose::PDB_file);
		core::pose::Pose pose = start_pose;



		core::conformation::symmetry::SymmData symmdata1(  pose.size(),  pose.num_jump() );
		std::string symm_def1 = "core/conformation/symmetry/symm_def1.dat";
		symmdata1.read_symmetry_data_from_file(symm_def1);
		core::pose::symmetry::make_symmetric_pose( pose, symmdata1 );

		TS_ASSERT( core::pose::symmetry::is_symmetric(pose));

		auto const & symm_conf (
			dynamic_cast<SymmetricConformation const & > ( pose.conformation() ) );
		SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

		AsymmetricUnitSelector master_selector = AsymmetricUnitSelector();
		utility::vector1< bool > subset = master_selector.apply(pose);

		for ( core::Size i = 1; i <= pose.size(); ++i ) {
			if ( symm_info->bb_is_independent(i) ) {
				TS_ASSERT( subset[i] == true);
			} else {
				TS_ASSERT( subset[i] == false);
			}
		}

	}
	void tearDown(){

	}







};
