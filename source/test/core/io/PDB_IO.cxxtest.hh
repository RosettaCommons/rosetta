// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/PDB_IO.cxxtest.hh
/// @brief  test suite for PDB reader/writer
/// @author Sergey Lyskov

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Package Headers
#include <core/io/pdb/pose_io.hh>
#include <core/import_pose/import_pose.hh>

// Project Headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>


using namespace core;

class PDB_IO : public CxxTest::TestSuite
{
	chemical::ResidueTypeSetCAP residue_set;

public:
	PDB_IO() {}

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-no_optH" );
		residue_set = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// ------------------------------------------ //
	// test that PDB input and output function properly
	void test_pdb_io() {
		pose::Pose pose;
		const std::string original_file_name("core/io/test_in.pdb");
		import_pose::pose_from_pdb(pose, original_file_name);

		// see if number of residues is correct
		TS_ASSERT_EQUALS( pose.total_residue(), 116u );

		// write pose to a new file...
		const std::string tmp_file_name("PDB_IO_cxxtest.pdb._tmp_");
		io::pdb::dump_pdb(pose, tmp_file_name);

		TS_ASSERT_FILE_EQ(original_file_name.c_str(), tmp_file_name.c_str());

		// read written file as new pose object
		pose::Pose P2;
		import_pose::pose_from_pdb(P2, tmp_file_name);

		// see if number of residues is correct
		TS_ASSERT_EQUALS( P2.total_residue(), 116u );


		// now test if residue information is the same
		bool should_exit = false;
		for(Size i=1; i<=pose.total_residue(); ++i) {
			conformation::Residue const & R1( pose.residue(i) );
			conformation::Residue const & R2( P2.residue(i) );

			// check if number of atoms is the same
			TS_ASSERT_EQUALS( R1.natoms(), R2.natoms());

			TS_ASSERT_EQUALS( R1.name3(), R2.name3());

			for(Size a=1; a<=R1.natoms(); a++) {
				conformation::Atom const & A1( R1.atom(a) );
				conformation::Atom const & A2( R2.atom(a) );

				TS_ASSERT_EQUALS( R1.atom_name(a), R2.atom_name(a));
				if( R1.atom_name(a) != R2.atom_name(a) ) should_exit=true;

				TS_ASSERT_EQUALS( A1.xyz().x(), A2.xyz().x());
				TS_ASSERT_EQUALS( A1.xyz().y(), A2.xyz().y());
				TS_ASSERT_EQUALS( A1.xyz().z(), A2.xyz().z());

				if( A1.xyz().x() != A2.xyz().x() ) should_exit=true;
				if( A1.xyz().y() != A2.xyz().y() ) should_exit=true;
				if( A1.xyz().z() != A2.xyz().z() ) should_exit=true;

				if( should_exit ) break;
			}
			if( should_exit ) break;
		}
	}

	void test_pdb_read_partial_residues() {
		pose::Pose pose;
		// This file has a leading fragment of an Arg residue that needs to be ignored.
		import_pose::pose_from_pdb(pose, "core/io/1ten.pdb");
		TS_ASSERT_EQUALS( pose.total_residue(), 89 );
		TS_TRACE(pose.annotated_sequence());
		TS_ASSERT_EQUALS( pose.annotated_sequence(),
				"L[LEU:NtermProteinFull]DAPSQIEVKDVTDTTALITWFKPLAEIDGIELTYGIKDVPGDRTTIDLTEDENQYSIGNLKPDTEYEVSLISRRGDMSSNPAKETFTT[THR:CtermProteinFull]" );
	}
};
