// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/mmcif.cxxtest.hh
/// @brief  test suite for basic mmcif reading/writing
/// @author Rocco Moretti

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Package Headers
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/mmcif/cif_writer.hh>
#include <core/io/mmcif/cif_reader.hh>
#include <core/io/StructFileReaderOptions.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pose_to_sfr/PoseToStructFileRepConverter.hh>

// Project Headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <utility/pointer/owning_ptr.hh>
#include <basic/Tracer.hh>

#include <cifparse/CifFile.h>
#include <cifparse/CifParserBase.h>
#include <cifparse/ISTable.h>

typedef utility::pointer::shared_ptr< CifFile > CifFileOP;
typedef utility::pointer::shared_ptr< CifParser > CifParserOP;

static THREAD_LOCAL basic::Tracer TR("core.io.mmcif_IO.cxxtest");

using namespace core;

class mmcif_IO : public CxxTest::TestSuite
{

public:
	mmcif_IO() {}

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-no_optH -obey_ENDMDL" );
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_mmcif_input() {
		// Of course cifparse isn't namespaced ...
		std::string diagnostics; //output from the parser about errors, etc.
		CifFileOP cifFile( new CifFile);
		{
			CifParserOP cifParser( new CifParser(cifFile.get()) );
			cifParser->Parse("core/io/1QYS.cif", diagnostics);
			if ( !diagnostics.empty() ) {
				TR.Error << "Diagnostics:" << std::endl;
				TR.Error << diagnostics << std::endl;
			}
			//cifParser no longer needed.
		}
		Block& block = cifFile->GetBlock(cifFile->GetFirstBlockName());
		ISTable& entry = block.GetTable("entry");
		std::string structure_id = entry(0, "id");
		TS_ASSERT_EQUALS(structure_id, "1QYS");

		ISTable& poly = block.GetTable("entity_poly");
		std::string seq = poly(0, "pdbx_seq_one_letter_code_can");
		TR << "Sequence " << seq << std::endl;
		TS_ASSERT_EQUALS(seq, "MGDIQVQVNIDDNGKNFDYTYTVTTESELQKVLNELMDYIKKQGAKRVRISITARTKKEAEKFAAILIKVFAELGYNDIN\nVTFDGDTVTVEGQLEGGSLEHHHHHH");

		ISTable& atoms = block.GetTable("atom_site");
		TS_ASSERT_EQUALS(atoms.GetNumRows(), 692 );
		TR << "Number of atoms: " << atoms.GetNumRows() << std::endl;
		std::string atom_2 = atoms(1, "label_atom_id");
		TS_ASSERT_EQUALS(atom_2, "CA" );
		TR << "Atom #2 is: " << atom_2 << std::endl;
	}

	void test_mmcif_vs_pdb_input() {
		core::pose::PoseOP pdb_pose = core::import_pose::pose_from_file( "core/io/1QYS.pdb", false , core::import_pose::PDB_file);
		core::pose::PoseOP cif_pose = core::import_pose::pose_from_file( "core/io/1QYS.cif", false , core::import_pose::CIF_file);;

		TR << "Total residue: ";
		TR << "pdb " << pdb_pose->total_residue() << "; ";
		TR << "cif " << cif_pose->total_residue() << std::endl;

		pdb_pose->dump_pdb( "from_pdb.pdb" );
		cif_pose->dump_pdb( "from_cif.pdb" );

		TS_ASSERT_EQUALS( pdb_pose->total_residue(), cif_pose->total_residue() );

		// Proxy for RT equality
		for ( Size i = 1; i <= pdb_pose->total_residue(); ++i ) {
			TS_ASSERT_EQUALS( pdb_pose->residue( i ).name(), cif_pose->residue( i ).name() );
		}

		TS_ASSERT_EQUALS( pdb_pose->fold_tree(), cif_pose->fold_tree() );

		// Check PDB Info here
		TS_ASSERT_EQUALS( pdb_pose->pdb_info()->chain(5), cif_pose->pdb_info()->chain(5) );
		TS_ASSERT_EQUALS( pdb_pose->pdb_info()->number(10), cif_pose->pdb_info()->number(10) );
	}

	void test_mmcif_output() {
		core::pose::PoseOP pdb_pose = core::import_pose::pose_from_file( "core/io/1QYS.pdb", false , core::import_pose::PDB_file);
		core::pose::PoseOP cif_pose = core::import_pose::pose_from_file( "core/io/1QYS.cif", false , core::import_pose::CIF_file);;

		TR << "Writing to cif: " << std::endl;
		core::io::StructFileReaderOptions sfro;
		core::io::StructFileRepOP pdb_sfr;
		core::io::pose_to_sfr::PoseToStructFileRepConverter pose_to_sfr;
		pose_to_sfr.init_from_pose( *pdb_pose );
		pose_to_sfr.sfr();
		pdb_sfr = pose_to_sfr.sfr();
		pose_to_sfr.new_sfr();
		pose_to_sfr.init_from_pose( *cif_pose );
		core::io::StructFileRepOP cif_sfr =  pose_to_sfr.sfr();

		core::io::mmcif::dump_cif( "from_pdb.cif", pdb_sfr, sfro );
		core::io::mmcif::dump_cif( "from_cif.cif", cif_sfr, sfro );

		TS_ASSERT_EQUALS( pdb_pose->fold_tree(), cif_pose->fold_tree() );

		// Check pose dump cif.
		TS_ASSERT_THROWS_NOTHING( pdb_pose->dump_cif( "from_pdb_pose_dump.cif") );

		// Check PDB Info here
		TS_ASSERT_EQUALS( pdb_pose->pdb_info()->chain(5), cif_pose->pdb_info()->chain(5) );
		TS_ASSERT_EQUALS( pdb_pose->pdb_info()->number(10), cif_pose->pdb_info()->number(10) );
	}

	void test_multimodel_cif() {
		core::pose::PoseOP pdb_pose = core::import_pose::pose_from_file( "core/io/1L8C.cif", false , core::import_pose::CIF_file);;
		TR << "Reading 1L8C from cif. " << std::endl;
		TS_ASSERT_EQUALS( pdb_pose->total_residue(), 149 );
	}

};
