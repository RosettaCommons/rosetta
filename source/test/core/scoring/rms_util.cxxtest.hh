// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

// C/C++ headers
#include <map>

// Project headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/id/SequenceMapping.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/util.hh>
#include <utility/vector1.hh>

using core::Real;
using core::Size;
using core::pose::Pose;
using core::pose::PoseOP;
using std::map;

class RmsUtilTest : public CxxTest::TestSuite {
public:
	void setUp() {
		// AMW: no longer have to separately load D amino acids
		core_init_with_additional_options( //Load the D-amino acids for this test.
			"-no_optH "//-extra_res_fa d-caa/DALA.params d-caa/DASP.params d-caa/DGLU.params d-caa/DPHE.params d-caa/DHIS.params d-caa/DHIS_D.params d-caa/DILE.params d-caa/DLYS.params d-caa/DLEU.params d-caa/DMET.params d-caa/DASN.params d-caa/DPRO.params d-caa/DGLN.params d-caa/DARG.params d-caa/DSER.params d-caa/DTHR.params d-caa/DVAL.params d-caa/DTRP.params d-caa/DTYR.params"
		);
		// Residue definitions can't be supplied on the command line b/c
		// the ResidueTypeSet is already initialized.
		using namespace core::chemical;
		utility::vector1< std::string > params_files;
		ResidueTypeSetCOP const_residue_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		ResidueTypeSet & residue_set = const_cast< ResidueTypeSet & >(*const_residue_set);
		if ( !residue_set.has_name("QC1") ) params_files.push_back("core/scoring/1pqc.params");
		residue_set.read_files_for_custom_residue_types(params_files);
	}

	void test_ligand_rms() {
		test::UTracer UT("core/scoring/rms_util.u");

		using namespace core;
		using namespace core::scoring;
		pose::Pose pose1, pose2;
		core::import_pose::pose_from_file( pose1, "core/scoring/1pqc_0001.pdb" , core::import_pose::PDB_file);
		core::import_pose::pose_from_file( pose2, "core/scoring/1pqc_0003.pdb" , core::import_pose::PDB_file);

		UT << "Ligand RMS, heavy atoms only:\n";
		UT << "    no superposition,   no automorphisms: " << rmsd_no_super(pose1, pose2, is_ligand_heavyatom) << "\n";
		UT << "  with superposition,   no automorphisms: " << rmsd_with_super(pose1, pose2, is_ligand_heavyatom) << "\n";
		UT << "    no superposition, with automorphisms: " << automorphic_rmsd(pose1.residue(1), pose2.residue(1), false /*don't superimpose*/) << "\n";
		UT << "  with superposition, with automorphisms: " << automorphic_rmsd(pose1.residue(1), pose2.residue(1), true /*superimpose*/) << "\n";
	}

	void test_gdtsc() {
		PoseOP ref = core::import_pose::pose_from_file("core/scoring/4dO8B.pdb", core::import_pose::PDB_file);
		PoseOP mod = core::import_pose::pose_from_file("core/scoring/model-5.pdb", core::import_pose::PDB_file);

		map<Size, Size> all_residues;
		for ( Size i = 1; i <= ref->size(); ++i ) {
			all_residues[i] = i;
		}

		TS_ASSERT_DELTA(0.745, core::scoring::gdtsc(*ref, *mod, all_residues), 0.01);
	}

	void test_gdtsc_self() {
		PoseOP mod = core::import_pose::pose_from_file("core/scoring/model-5.pdb", core::import_pose::PDB_file);

		map<Size, Size> all_residues;
		for ( Size i = 1; i <= mod->size(); ++i ) {
			all_residues[i] = i;
		}

		TS_ASSERT_DELTA(1.0, core::scoring::gdtsc(*mod, *mod, all_residues), 0.0001);
	}

	void test_superimpose_self() {
		Pose pose = *core::import_pose::pose_from_file("core/scoring/2GB3.pdb", core::import_pose::PDB_file);
		Real rmsd = core::scoring::CA_rmsd(pose, pose);
		TS_ASSERT_DELTA(0, rmsd, 0.0001);
	}

	void test_superimpose_same() {
		Pose pose1 = *core::import_pose::pose_from_file("core/scoring/2GB3.pdb", core::import_pose::PDB_file);
		Pose pose2 = *core::import_pose::pose_from_file("core/scoring/2GB3.pdb", core::import_pose::PDB_file);
		Real rmsd = core::scoring::CA_rmsd(pose1, pose2);
		TS_ASSERT_DELTA(0, rmsd, 0.0001);
	}

	void test_ca_rmsd_with_mapping() {
		Pose pose1 = *core::import_pose::pose_from_file("core/scoring/2GB3.pdb", core::import_pose::PDB_file);
		Pose pose2 = *core::import_pose::pose_from_file("core/scoring/2GB3.pdb", core::import_pose::PDB_file);

		std::map<Size, Size> residues;
		for ( Size i = 1; i <= pose1.size(); ++i ) {
			residues[i] = i;
		}

		TS_ASSERT_DELTA(0, core::scoring::CA_rmsd(pose1, pose2, residues), 0.0001);
	}

	/// @brief Test rmsd calculation as applied to a mixed D/L helical bundle.
	/// @details This uses a parametrically-generated helical bundle as the native pose,
	/// with aligned and unaligned helical bundles (CA RMSD=9.921 A) that are supposed to
	/// be aligned to the native.  This test asks whether the aligned and unaligned
	/// poses yield the same RMSD.
	/// @author Vikram K. Mulligan (vmullig@uw.edu), Baker Laboratory
	void test_DLbundle_rmsd_aligned_unaligned()
	{
		PoseOP pose_native = core::import_pose::pose_from_file( "core/scoring/native_for_rmstest.pdb", false , core::import_pose::PDB_file);
		PoseOP pose_unaligned = core::import_pose::pose_from_file( "core/scoring/unaligned_for_rmstest.pdb", false , core::import_pose::PDB_file);
		PoseOP pose_aligned = core::import_pose::pose_from_file( "core/scoring/aligned_for_rmstest.pdb", false , core::import_pose::PDB_file);

		Real rmsd_unaligned = core::scoring::CA_rmsd( *pose_unaligned, *pose_native );
		Real rmsd_aligned = core::scoring::CA_rmsd( *pose_aligned, *pose_native );

		TS_ASSERT_DELTA( rmsd_unaligned, rmsd_aligned, 0.0001 );

		return;
	} //test_DLbundle_rmsd_aligned_unaligned()

};
