// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/PDB_IO.cxxtest.hh
/// @brief  test suite for PDB reader/writer
/// @author Sergey Lyskov

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Package Headers
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/StructFileRepOptions.fwd.hh>
#include <core/import_pose/import_pose.hh>

// Project Headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>
#include <fstream>

static basic::Tracer TR("core.io.PDB_IO.cxxtest");

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
		import_pose::pose_from_file(pose, original_file_name, core::import_pose::PDB_file);

		// see if number of residues is correct
		TS_ASSERT_EQUALS( pose.size(), 116u );

		// write pose to a new file...
		const std::string tmp_file_name("PDB_IO_cxxtest.pdb._tmp_");
		io::pdb::dump_pdb(pose, tmp_file_name);

		TS_ASSERT_FILE_EQ(original_file_name.c_str(), tmp_file_name.c_str());

		// read written file as new pose object
		pose::Pose P2;
		import_pose::pose_from_file(P2, tmp_file_name, core::import_pose::PDB_file);

		// see if number of residues is correct
		TS_ASSERT_EQUALS( P2.size(), 116u );


		// now test if residue information is the same
		bool should_exit = false;
		for ( Size i=1; i<=pose.size(); ++i ) {
			conformation::Residue const & R1( pose.residue(i) );
			conformation::Residue const & R2( P2.residue(i) );

			// check if number of atoms is the same
			TS_ASSERT_EQUALS( R1.natoms(), R2.natoms());

			TS_ASSERT_EQUALS( R1.name3(), R2.name3());

			for ( Size a=1; a<=R1.natoms(); a++ ) {
				conformation::Atom const & A1( R1.atom(a) );
				conformation::Atom const & A2( R2.atom(a) );

				TS_ASSERT_EQUALS( R1.atom_name(a), R2.atom_name(a));
				if ( R1.atom_name(a) != R2.atom_name(a) ) should_exit=true;

				TS_ASSERT_EQUALS( A1.xyz().x(), A2.xyz().x());
				TS_ASSERT_EQUALS( A1.xyz().y(), A2.xyz().y());
				TS_ASSERT_EQUALS( A1.xyz().z(), A2.xyz().z());

				if ( A1.xyz().x() != A2.xyz().x() ) should_exit=true;
				if ( A1.xyz().y() != A2.xyz().y() ) should_exit=true;
				if ( A1.xyz().z() != A2.xyz().z() ) should_exit=true;

				if ( should_exit ) break;
			}
			if ( should_exit ) break;
		}
	}

	// make sure that SEQRES records are written when the write_seqres_records option is specified
	void test_pdb_io_seqres() {
		pose::Pose pose;
		const std::string original_file_name("core/io/test_in.pdb");
		import_pose::pose_from_file(pose, original_file_name, core::import_pose::PDB_file);

		core::io::StructFileRepOptionsOP options( new core::io::StructFileRepOptions );
		options->set_write_seqres_records(true);

		std::stringstream output_data;
		io::pdb::dump_pdb(pose, output_data, core::io::StructFileRepOptionsCOP(options));

		utility::vector1<std::string> expected_seqres {
			"SEQRES   1 A  116  ASP ALA ILE THR ILE HIS SER ILE LEU ASP TRP ILE GLU          ",
			"SEQRES   2 A  116  ASP ASN LEU GLU SER PRO LEU SER LEU GLU LYS VAL SER          ",
			"SEQRES   3 A  116  GLU ARG SER GLY TYR SER LYS TRP HIS LEU GLN ARG MET          ",
			"SEQRES   4 A  116  PHE LYS LYS GLU THR GLY HIS SER LEU GLY GLN TYR ILE          ",
			"SEQRES   5 A  116  ARG SER ARG LYS MET THR GLU ILE ALA GLN LYS LEU LYS          ",
			"SEQRES   6 A  116  GLU SER ASN GLU PRO ILE LEU TYR LEU ALA GLU ARG TYR          ",
			"SEQRES   7 A  116  GLY PHE GLU SER GLN GLN THR LEU THR ARG THR PHE LYS          ",
			"SEQRES   8 A  116  ASN TYR PHE ASP VAL PRO PRO HIS LYS TYR ARG MET THR          ",
			"SEQRES   9 A  116  ASN MET GLN GLY GLU SER ARG PHE LEU HIS PRO LEU              "};

		// std::ifstream infile(tmp_file_name);
		std::string file_line;
		utility::vector1<std::string> actual_seqres;
		while ( std::getline(output_data, file_line) ) {
			if ( file_line.substr(0,6) == "SEQRES" ) {
				actual_seqres.push_back(file_line);
			}
		}

		TS_ASSERT_EQUALS(actual_seqres, expected_seqres);
	}

	void test_pdb_read_partial_residues() {
		pose::Pose pose;
		// This file has a leading fragment of an Arg residue that needs to be ignored.
		import_pose::pose_from_file(pose, "core/io/1ten.pdb", core::import_pose::PDB_file);
		TS_ASSERT_EQUALS( pose.size(), 89 );
		TR << pose.annotated_sequence() << std::endl;
		TS_ASSERT_EQUALS( pose.annotated_sequence(),
			"L[LEU:NtermProteinFull]DAPSQIEVKDVTDTTALITWFKPLAEIDGIELTYGIKDVPGDRTTIDLTEDENQYSIGNLKPDTEYEVSLISRRGDMSSNPAKETFTT[THR:CtermProteinFull]" );
	}
};
