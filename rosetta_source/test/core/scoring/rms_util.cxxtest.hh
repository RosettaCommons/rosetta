// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>

using core::Real;
using core::pose::Pose;

class RmsUtilTest : public CxxTest::TestSuite {
 public:
	void setUp() {
		core_init_with_additional_options(
			"-no_optH"// -extra_res_fa core/scoring/1pqc.params"
		);
		// Residue definitions can't be supplied on the command line b/c
		// the ResidueTypeSet is already initialized.
		using namespace core::chemical;
		utility::vector1< std::string > params_files;
		ResidueTypeSetCAP const_residue_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		ResidueTypeSet & residue_set = const_cast< ResidueTypeSet & >(*const_residue_set);
		if(!residue_set.has_name("QC1")) params_files.push_back("core/scoring/1pqc.params");
		residue_set.read_files(params_files,
			ChemicalManager::get_instance()->atom_type_set( FA_STANDARD ),
			ChemicalManager::get_instance()->element_set( FA_STANDARD ),
			ChemicalManager::get_instance()->mm_atom_type_set( FA_STANDARD ),
			ChemicalManager::get_instance()->orbital_type_set(FA_STANDARD));//,
	}

	void test_ligand_rms() {
		test::UTracer UT("core/scoring/rms_util.u");

		using namespace core;
		using namespace core::scoring;
		pose::Pose pose1, pose2;
		core::import_pose::pose_from_pdb( pose1, "core/scoring/1pqc_0001.pdb" );
		core::import_pose::pose_from_pdb( pose2, "core/scoring/1pqc_0003.pdb" );

		UT << "Ligand RMS, heavy atoms only:\n";
		UT << "    no superposition,   no automorphisms: " << rmsd_no_super(pose1, pose2, is_ligand_heavyatom) << "\n";
		UT << "  with superposition,   no automorphisms: " << rmsd_with_super(pose1, pose2, is_ligand_heavyatom) << "\n";
		UT << "    no superposition, with automorphisms: " << automorphic_rmsd(pose1.residue(1), pose2.residue(1), false /*don't superimpose*/) << "\n";
		UT << "  with superposition, with automorphisms: " << automorphic_rmsd(pose1.residue(1), pose2.residue(1), true /*superimpose*/) << "\n";
	}

  void test_superimpose_self() {
    Pose pose = *core::import_pose::pose_from_pdb("core/scoring/2GB3.pdb");
    Real rmsd = core::scoring::CA_rmsd(pose, pose);
    TS_ASSERT_DELTA(0, rmsd, 0.0001);
  }

  void test_superimpose_same() {
    Pose pose1 = *core::import_pose::pose_from_pdb("core/scoring/2GB3.pdb");
    Pose pose2 = *core::import_pose::pose_from_pdb("core/scoring/2GB3.pdb");
    Real rmsd = core::scoring::CA_rmsd(pose1, pose2);
    TS_ASSERT_DELTA(0, rmsd, 0.0001);
  }
};
