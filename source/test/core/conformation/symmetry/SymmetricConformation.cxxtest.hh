// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/optimization/Minimizer.cxxtest.hh
/// @brief  test suite for Minimizer
/// @author Phil Bradley
/// @author Sergey Lyskov

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>

// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

// AUTO-REMOVED #include <core/import_pose/import_pose.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>

#include <core/scoring/rms_util.hh>

#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmData.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.conformation.symmetry.SymmetricConformation.cxxtest");

using namespace core;

class SymmetricConformationTests : public CxxTest::TestSuite
{
	chemical::ResidueTypeSetCAP residue_set;

public:
	SymmetricConformationTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();
		residue_set = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );
	}

	// Shared finalization goes here.
	void tearDown() {
	}


///////////////////////////////////////////////////////////////////////////////
// ------------------------------------------ //
/// @brief simple test minimization
void test_symmetric_conformation()
{
	double rms_threshold = 1e-2;

	pose::Pose start_pose;
	import_pose::pose_from_pdb( start_pose, "core/scoring/symmetry/test_in.pdb" );
	pose::Pose pose = start_pose;

	// test different symmetry data input files
	core::conformation::symmetry::SymmData symmdata1(  pose.n_residue(),  pose.num_jump() );
	std::string symm_def1 = "core/conformation/symmetry/symm_def1.dat";
	symmdata1.read_symmetry_data_from_file(symm_def1);
	core::pose::symmetry::make_symmetric_pose( pose, symmdata1 );
	pose::Pose pose1_ref;
	import_pose::pose_from_pdb( pose1_ref, "core/scoring/symmetry/symm_test.pdb" );
	Real rms_to_restored = scoring::CA_rmsd( pose, pose1_ref );
	TR << "RMS difference after symmetry reconstruction of test1 trimer from monomer : " << rms_to_restored << std::endl;
	  TS_ASSERT( rms_to_restored < rms_threshold );

	import_pose::pose_from_pdb( pose, "core/conformation/symmetry/symm_test_in2.pdb" );;
	core::conformation::symmetry::SymmData symmdata2(  pose.n_residue(),  pose.num_jump() );
	std::string symm_def2 = "core/conformation/symmetry/symm_def2.dat";
	symmdata2.read_symmetry_data_from_file(symm_def2);
	core::pose::symmetry::make_symmetric_pose( pose, symmdata2 );
	pose::Pose pose2_ref;
	import_pose::pose_from_pdb( pose2_ref, "core/conformation/symmetry/symm_test2.pdb" );
	rms_to_restored = scoring::CA_rmsd( pose, pose2_ref );
	TR << "RMS difference after symmetry reconstruction of test2 D6 from monomer : " << rms_to_restored << std::endl;
	  TS_ASSERT( rms_to_restored < rms_threshold );

};


};
