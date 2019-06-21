// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/mmtf.cxxtest.hh
/// @brief  test suite for basic mmtf reading/writing
/// @author Danny Farrell

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Package Headers
#include <core/chemical/AA.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <src/core/io/mmtf/mmtf_reader.hh>
#include <src/core/io/mmtf/mmtf_writer.hh>
#include <src/core/io/mmcif/cif_reader.hh>
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

#include <mmtf.hpp>
#include <cifparse/CifFile.h>
#include <cifparse/CifParserBase.h>

static basic::Tracer TR("core.io.mmtf_IO.cxxtest");

using namespace core;

class mmtf_IO : public CxxTest::TestSuite
{
private:
	core::pose::PoseOP pdb_pose;
	utility::vector1< std::string > mmtf_test_files;

	// Helper function three letter code to one letter code (doesnt fail, returns X)
	char tlc_2_olc( std::string const & name )
	{
		TR << "tlc 2 olc: " << name << std::endl;
		if ( core::chemical::is_aa_name_unknown( name ) ) return 'X';
		return core::chemical::oneletter_code_from_aa( core::chemical::aa_from_one_or_three( name ) );
	}

	std::string atoms_to_sequence( core::io::ChainAtoms const & all_atoms ) {
		std::string sequence("");
		for ( auto const & atm : all_atoms ) sequence += tlc_2_olc( atm.resName );
		return sequence;
	}

public:
	mmtf_IO() {}

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options("-ignore_zero_occupancy false -ignore_unrecognized_res -packing::pack_missing_sidechains false");
		pdb_pose = core::import_pose::pose_from_file( "core/io/1QYS.pdb", false , core::import_pose::PDB_file);
		// "173D", "1BNA", // missing group name
		// "1IGT", // has BMA hetatm
		mmtf_test_files.push_back("1MSH");
		mmtf_test_files.push_back("1SKM");
		mmtf_test_files.push_back("4CK4");
		mmtf_test_files.push_back("4P3R");
		mmtf_test_files.push_back("5EMG"); // "empty-mmtfVersion99999999.mmtf",
		mmtf_test_files.push_back("1AA6");
		mmtf_test_files.push_back("1CAG");
		mmtf_test_files.push_back("1L2Q");
		mmtf_test_files.push_back("1O2F");
		mmtf_test_files.push_back("3NJW");
		mmtf_test_files.push_back("4CUP");
		// "4V5A", // too large for pdb file comparison
		mmtf_test_files.push_back("5ESW"); // "empty-numChains1.mmtf",
		mmtf_test_files.push_back("1AUY");
		mmtf_test_files.push_back("1HTQ");
		mmtf_test_files.push_back("1LPV");
		// mmtf_test_files.push_back("1R9V"); // failed to set abase2 for acceptor atom, it has no nbrs! (pdb fails)
		// mmtf_test_files.push_back("3ZYB"); // mmcif mmtf difference
		mmtf_test_files.push_back("4OPJ");
		mmtf_test_files.push_back("4Y60");
		// "empty-all0.mmtf", "empty-numModels1.mmtf"
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	// Read a pdb file, write it to a mmtf file, read that mmtf file, and then write it again, and then compare the two
	// mmtf files to make sure they're the same.
	void test_mmtf_roundtrip() {
		core::Size initial_size = pdb_pose->size();

		std::string mmtf_file_1 = "first_out.mmtf";
		pdb_pose->dump_mmtf(mmtf_file_1);

		core::pose::PoseOP mmtf_pose_1 = pose_from_file(mmtf_file_1, false,
			core::import_pose::FileType::MMTF_file);
		TS_ASSERT_EQUALS( initial_size, mmtf_pose_1->size() );

		std::string mmtf_file_2 = "second_out.mmtf";
		mmtf_pose_1->dump_mmtf(mmtf_file_2);

		core::pose::PoseOP mmtf_pose_2 = pose_from_file(mmtf_file_2, false,
			core::import_pose::FileType::MMTF_file);
		TS_ASSERT_EQUALS( initial_size, mmtf_pose_2->size() );

		mmtf::StructureData sd1;
		mmtf::decodeFromFile(sd1, mmtf_file_1);
		mmtf::StructureData sd2;
		mmtf::decodeFromFile(sd2, mmtf_file_2);
		TS_ASSERT_EQUALS( true, sd1 == sd2 );
	}

	/** we cannot check sd1 == sd2 because rosetta will add atoms
	*  so instead we simply compare a cif import to an mmtf import
	*  for each cif file in the MMTF test suite and see if they have
	*  the same # of residues, and that sequence is the same.
	*/
	void test_mmtf_testsuite() {
		for ( auto const & test_file : mmtf_test_files ) {
			TR << "Working on mmtf test file: " << test_file << std::endl;
			std::string mmtf_file = "core/io/mmtf/" + test_file + ".mmtf";
			std::string cif_file = "core/io/mmtf/" + test_file + ".cif";
			mmtf::StructureData sd1;
			mmtf::decodeFromFile(sd1, mmtf_file);
			TR << "Loaded files here" << std::endl;
			core::io::StructFileReaderOptions opts;

			core::io::StructFileRepOP mmtfsfr(core::io::mmtf::create_sfr_from_mmtf_filename( mmtf_file, opts ));
			TR << "Loaded mmtf file" << std::endl;

			std::string contents_of_file;
			utility::io::izstream file( cif_file );
			utility::slurp( file, contents_of_file );
			std::string diagnostics;
			CifFileOP cifFile( new CifFile );
			CifParserOP cifParser( new CifParser( cifFile.get() ) );
			cifParser->ParseString( contents_of_file, diagnostics );
			core::io::StructFileRepOP cifsfr( core::io::mmcif::create_sfr_from_cif_file_op( cifFile, opts ) );
			TR << "Loaded cif file" << std::endl;

			TS_ASSERT_EQUALS( cifsfr->chains().size(), mmtfsfr->chains().size() );
			for ( core::Size i=0; i < cifsfr->chains().size(); ++i ) {
				TS_ASSERT_EQUALS( cifsfr->chains()[i].size(), mmtfsfr->chains()[i].size() );
				utility::vector1< std::string > mmtfnames, cifnames;
				for ( auto const & x : mmtfsfr->chains()[i] ) mmtfnames.push_back(x.resName);
				for ( auto const & x : cifsfr->chains()[i] ) cifnames.push_back(x.resName);
				TS_ASSERT_EQUALS(mmtfnames, cifnames);
			}
		}
	}

	void test_aiPose_from_sfr() {
		core::io::StructFileRepOptionsOP options =  core::io::StructFileRepOptionsOP( new core::io::StructFileRepOptions );
		core::io::mmtf::set_mmtf_default_options( *options );
		core::io::pose_to_sfr::PoseToStructFileRepConverter converter = core::io::pose_to_sfr::PoseToStructFileRepConverter( *options );
		converter.init_from_pose( *pdb_pose );
		core::io::StructFileRepOP sfr =  converter.sfr();

		core::Size chain_num = 0, group_num = 1;
		auto ai_pose( core::io::mmtf::aiPose_from_sfr(*sfr) );

		TS_ASSERT_EQUALS(true, ai_pose[0][0][0].resName == pdb_pose->residue(1).name3() );

		for ( core::Size i=2; i<=pdb_pose->size(); ++i ) {
			char prev_chain_id( pdb_pose->pdb_info()->chain(i-1) );
			char chain_id( pdb_pose->pdb_info()->chain(i) );
			if ( prev_chain_id != chain_id ) {
				++chain_num;
				group_num = 0;
			}
			auto ai = ai_pose[chain_num][group_num][0];
			TS_ASSERT_EQUALS(true, ai.resName == pdb_pose->residue(i).name3() );
			TS_ASSERT_EQUALS(true, ai_pose[chain_num][group_num].size() == pdb_pose->residue(i).natoms());
			++group_num;
		}
	}
};
