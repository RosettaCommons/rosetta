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
		core_init_with_additional_options( "-no_optH -INTEGRATION_TEST" );
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

	void test_pdb_read_d_aa() {
		using namespace core::chemical;

		pose::Pose pose;
		// Input pose has the 20 Rosetta three letter D-AA codes, followed by the 6 D-AA entries using the wwPDB code which is different from the Rosetta code,
		// followed by the 6 wwPDB entries (as separate ligand chains) which overlap with the Rosetta D-AA code.
		// The first 20+6 entries should get loaded as Rosetta D-AA types.
		// The final 6 entries don't necessarily need to get loaded as CCD components, but they *shouldn't* be loaded as the D-AA which just happens to have the same code.
		// (These are currently disabled because of PDB loading issues.)
		import_pose::pose_from_file(pose, "core/io/daa.pdb", core::import_pose::PDB_file);
		TS_ASSERT_EQUALS( pose.size(), 26 );
		//TS_ASSERT_EQUALS( pose.size(), 32 );
		TR << pose.annotated_sequence() << std::endl;
		TS_ASSERT_EQUALS( pose.annotated_sequence(),
			"A[DALA:NtermProteinFull]C[DCYS]D[DASP]E[DGLU]F[DPHE]GH[DHIS]I[DILE]K[DLYS]L[DLEU]M[DMET]N[DASN]P[DPRO]Q[DGLN]R[DARG]S[DSER]T[DTHR]V[DVAL]W[DTRP]Y[DTYR]"
			"C[DCYS]E[DGLU]F[DPHE]M[DMET]D[DASP]S[DSER:CtermProteinFull]"
			//"Z[pdb_DCS]Z[pdb_DGU]Z[pdb_DPH]Z[pdb_DME]Z[pdb_DAS]Z[pdb_DSE:CtermProteinFull:NtermProteinFull]"
		);
		utility::vector1< std::string > three_letter_codes = { "DAL", "DCS", "DAS", "DGU", "DPH", "GLY", "DHI", "DIL", "DLY", "DLE", "DME", "DAN", "DPR", "DGN", "DAR", "DSE", "DTH", "DVA", "DTR", "DTY", // Rosetta specific codes
			"DCS", "DGU", "DPH", "DME", "DAS", "DSE" // Loaded as CCD codes, but in-pose as Rosetta ones.
			//"DCS", "DGU", "DPH", "DME", "DAS", "DSE" // The (non D-AA) CCD entries
			};

		for ( core::Size ii(1); ii <= pose.size(); ++ii ) {
			TS_ASSERT_EQUALS( three_letter_codes[ii], pose.residue(ii).name3() );
		}


		// Iterate over the Rosetta residue types, double checking their properties are as expected.
		utility::vector1< AA > const aa_enums = { aa_dal, aa_dcs, aa_das, aa_dgu, aa_dph, aa_gly, aa_dhi, aa_dil, aa_dly, aa_dle, aa_dme, aa_dan, aa_dpr, aa_dgn, aa_dar, aa_dse, aa_dth, aa_dva, aa_dtr, aa_dty,
			aa_dcs, aa_dgu, aa_dph, aa_dme, aa_das, aa_dse };
		utility::vector1< std::string > basenames = { "DALA", "DCYS", "DASP", "DGLU", "DPHE", "GLY", "DHIS", "DILE", "DLYS", "DLEU", "DMET", "DASN", "DPRO", "DGLN", "DARG", "DSER", "DTHR", "DVAL", "DTRP", "DTYR",
			"DCYS", "DGLU", "DPHE", "DMET", "DASP", "DSER" };
		utility::vector1< bool > const ispolar = { false, false, true, true, false, false, true, false, true, false, false, true, false, true, true, true, true, false, false, false,
			false, true, false, false, true, true };

		TS_ASSERT_EQUALS( aa_enums.size(), basenames.size() );
		TS_ASSERT_EQUALS( aa_enums.size(), ispolar.size() );

		for ( core::Size ii(1); ii <= aa_enums.size(); ++ii ) {
			TS_ASSERT_EQUALS( aa_enums[ii], pose.residue_type(ii).aa() );
			TS_ASSERT_EQUALS( basenames[ii], pose.residue_type(ii).base_name() );
			TS_ASSERT_EQUALS( bool(ispolar[ii]), bool(pose.residue_type(ii).is_polar()) ); // Odd compile error if I don't cast

			TS_ASSERT( pose.residue_type(ii).is_alpha_aa() );
			if ( pose.residue_type(ii).base_name() != "GLY" ) {
				TS_ASSERT( pose.residue_type(ii).is_d_aa() );
			}
		}
	}
};
