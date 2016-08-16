// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/PDBInfo.cxxtest.hh
/// @brief  test suite for core::pose::PDBInfo
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// Test headers
#include <cxxtest/TestSuite.h>

// Package headers
#include <test/core/init_util.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/Pose.hh>

// utility headers
#include <utility/vector1.hh>

//Auto Headers
#include <core/id/AtomID_Mask.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/AtomPointer.hh>
#include <core/kinematics/DomainMap.hh>
#include <utility/io/all.fwd.hh>
#include <utility/io/icstream.fwd.hh>
#include <utility/io/ocstream.fwd.hh>
#include <numeric/BodyPosition.fwd.hh>
#include <numeric/Quaternion.fwd.hh>
#include <numeric/all.fwd.hh>


// --------------- Test Class --------------- //

class PDBInfoTests : public CxxTest::TestSuite {

public:

	// typedefs
	typedef core::Size Size;
	typedef core::conformation::Conformation Conformation;
	typedef core::conformation::ConformationOP ConformationOP;
	typedef core::conformation::ResidueOP ResidueOP;
	typedef core::pose::PDBInfo PDBInfo;
	typedef core::pose::PDBPoseMap PDBPoseMap;
	typedef core::pose::Pose Pose;

	// shared data
	Pose pose;
	ResidueOP ala_rsd;

	// shared initialization
	void setUp() {
		using namespace core::chemical;
		using namespace core::conformation;

		core_init();
		core::import_pose::pose_from_file( pose, "core/pose/pdbinfo_test_in.pdb" , core::import_pose::PDB_file);

		ResidueTypeSetCOP residue_set
			( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
		ala_rsd = ResidueFactory::create_residue( residue_set->name_map( "ALA" ) );
	}

	// shared finalization
	void tearDown() {
	}

	// --------------- Test Cases --------------- //

	/// @brief test basic functionality of PDBInfo from input
	void test_PDBInfo_from_input_pdb() {
		// chain
		TS_ASSERT_EQUALS( pose.pdb_info()->chain( 1 ), 'L' );
		TS_ASSERT_EQUALS( pose.pdb_info()->chain( 2 ), 'H' );
		TS_ASSERT_EQUALS( pose.pdb_info()->chain( 3 ), 'H' );
		TS_ASSERT_EQUALS( pose.pdb_info()->chain( 4 ), 'P' );

		// residue number
		TS_ASSERT_EQUALS( pose.pdb_info()->number( 1 ), 1 );
		TS_ASSERT_EQUALS( pose.pdb_info()->number( 2 ), 100 );
		TS_ASSERT_EQUALS( pose.pdb_info()->number( 3 ), 100 );
		TS_ASSERT_EQUALS( pose.pdb_info()->number( 4 ), 671 );

		// insertion code
		TS_ASSERT_EQUALS( pose.pdb_info()->icode( 1 ), ' ' );
		TS_ASSERT_EQUALS( pose.pdb_info()->icode( 2 ), 'C' );
		TS_ASSERT_EQUALS( pose.pdb_info()->icode( 3 ), 'D' );
		TS_ASSERT_EQUALS( pose.pdb_info()->icode( 4 ), ' ' );
	}

	/// @brief test basic functionality of internally maintained PDBPoseMap
	///  of PDBInfo from input
	void test_PDBInfo_from_input_pdb_PDBPoseMap() {
		TS_ASSERT_EQUALS( pose.pdb_info()->pdb2pose().size(), 4 );
		TS_ASSERT_EQUALS( pose.pdb_info()->pdb2pose( 'L', 1, ' ' ), 1 );
		TS_ASSERT_EQUALS( pose.pdb_info()->pdb2pose( 'H', 100, 'C' ), 2 );
		TS_ASSERT_EQUALS( pose.pdb_info()->pdb2pose( 'H', 100, 'D' ), 3 );
		TS_ASSERT_EQUALS( pose.pdb_info()->pdb2pose( 'P', 671, ' ' ), 4 );
	}

	/// @brief test single residue chain mutator
	void test_PDBInfo_single_residue_chain_mutator() {
		// make a copy of the input pdb_info
		PDBInfo info = *pose.pdb_info();

		// invoke mutator
		info.chain( 1, 'A' );
		TS_ASSERT_EQUALS( info.chain( 1 ), 'A' );
		TS_ASSERT_EQUALS( info.pdb2pose( 'A', 1, ' ' ), 1 );
		TS_ASSERT_EQUALS( info.pdb2pose( 'L', 1, ' ' ), 0 ); // 0 == not found
	}

	/// @brief test single residue number mutator
	void test_PDBInfo_single_residue_number_mutator() {
		// make a copy of the input pdb_info
		PDBInfo info = *pose.pdb_info();

		// invoke mutator
		info.number( 2, 999 );
		TS_ASSERT_EQUALS( info.number( 2 ), 999 );
		TS_ASSERT_EQUALS( info.pdb2pose( 'H', 999, 'C' ), 2 );
		TS_ASSERT_EQUALS( info.pdb2pose( 'H', 100, 'C' ), 0 ); // 0 == not found
	}

	/// @brief test single residue insertion code mutator
	void test_PDBInfo_single_residue_icode_mutator() {
		// make a copy of the input pdb_info
		PDBInfo info = *pose.pdb_info();

		// invoke mutator
		info.icode( 3, 'Z' );
		TS_ASSERT_EQUALS( info.icode( 3 ), 'Z' );
		TS_ASSERT_EQUALS( info.pdb2pose( 'H', 100, 'Z' ), 3 );
		TS_ASSERT_EQUALS( info.pdb2pose( 'H', 100, 'D' ), 0 ); // 0 == not found
	}

	/// @brief test single residue resinfo mutator
	void test_PDBInfo_single_residue_resinfo_mutator() {
		// make a copy of the input pdb_info
		PDBInfo info = *pose.pdb_info();

		// invoke mutator
		info.set_resinfo( 4, 'X', -1, 'Y' );
		TS_ASSERT_EQUALS( info.chain( 4 ), 'X' );
		TS_ASSERT_EQUALS( info.number( 4 ), -1 );
		TS_ASSERT_EQUALS( info.icode( 4 ), 'Y' );
		TS_ASSERT_EQUALS( info.pdb2pose( 'X', -1, 'Y' ), 4 );
		TS_ASSERT_EQUALS( info.pdb2pose( 'P', 671, ' ' ), 0 ); // 0 == not found
	}

	/// @brief test en masse residue mutators
	void test_PDBInfo_en_masse_residue_mutators() {
		// make a copy of the input pdb_info
		PDBInfo info = *pose.pdb_info();

		// make fake chain/resid/icode arrays
		utility::vector1< char > chain;
		chain.push_back( 'A' );
		chain.push_back( 'B' );
		chain.push_back( 'C' );
		chain.push_back( 'D' );
		utility::vector1< int > resid;
		resid.push_back( -2 );
		resid.push_back( -1 );
		resid.push_back( 0 );
		resid.push_back( 1 );
		utility::vector1< char > icode;
		icode.push_back( 'W' );
		icode.push_back( 'X' );
		icode.push_back( 'Y' );
		icode.push_back( 'Z' );

		// invoke mutators
		info.set_chains( chain );
		info.set_numbering( resid );
		info.set_icodes( icode );

		TS_ASSERT_EQUALS( info.chain( 1 ), 'A' );
		TS_ASSERT_EQUALS( info.number( 2 ), -1 );
		TS_ASSERT_EQUALS( info.icode( 4 ), 'Z' );
		TS_ASSERT_EQUALS( info.pdb2pose( 'A', -2, 'W' ), 1 );
		TS_ASSERT_EQUALS( info.pdb2pose( 'H', 100, 'C' ), 0 ); // 0 == not found
	}

	/// @brief test residue+atom record resize and availability
	void test_PDBInfo_resize() {
		PDBInfo info;
		info.resize_residue_records( 4 );
		info.resize_atom_records( pose );
		info.occupancy( 4, 1, 0.5 );

		TS_ASSERT_EQUALS( info.natoms( 3 ), 8 ); // 7 atoms in gly + 1 OXT atom
		TS_ASSERT_EQUALS( info.occupancy( 4, 1 ), 0.5 );
	}

	/// @brief test append residue
	void test_PDBInfo_append_residue() {
		PDBInfo info = *pose.pdb_info();

		info.append_res( 3, 7, 2 ); // 2 residues after resid 3 with 7 atoms apiece

		TS_ASSERT_EQUALS( info.nres(), 6 );
		TS_ASSERT_EQUALS( info.natoms( 4 ), 7 );
		TS_ASSERT_EQUALS( info.natoms( 5 ), 7 );
		TS_ASSERT_EQUALS( info.chain( 5 ), PDBInfo::empty_record() );
		TS_ASSERT_EQUALS( info.chain( 3 ), 'H' );
		TS_ASSERT_EQUALS( info.chain( 6 ), 'P' );
	}

	/// @brief test prepend residue
	void test_PDBInfo_prepend_residue() {
		PDBInfo info = *pose.pdb_info();

		info.prepend_res( 3, 7, 2 ); // 2 residues before resid 3 with 7 atoms apiece

		TS_ASSERT_EQUALS( info.nres(), 6 );
		TS_ASSERT_EQUALS( info.natoms( 3 ), 7 );
		TS_ASSERT_EQUALS( info.natoms( 4 ), 7 );
		TS_ASSERT_EQUALS( info.chain( 3 ), PDBInfo::empty_record() );
		TS_ASSERT_EQUALS( info.chain( 2 ), 'H' );
		TS_ASSERT_EQUALS( info.chain( 5 ), 'H' );
	}

	/// @brief test delete residue
	void test_PDBInfo_delete_residue() {
		PDBInfo info = *pose.pdb_info();

		info.delete_res( 2 );
		TS_ASSERT_EQUALS( info.nres(), 3 );
		TS_ASSERT_EQUALS( info.natoms( 2 ), 8 ); // 7 atoms in gly + 1 OXT atom
		TS_ASSERT_EQUALS( info.chain( 3 ), 'P' );
	}

	/// @brief test observer attach to/detach from Conformation
	void test_PDBInfo_observe_attach_detach() {
		Conformation & conf = pose.conformation();
		PDBInfo info = *pose.pdb_info();

		info.attach_to( conf );
		TS_ASSERT_EQUALS( info.is_observing().lock().get(), &conf );
		info.detach_from();
	}

	/// @brief test Conformation observer behavior on append
	void test_PDBInfo_observe_append() {
		Conformation & conf = pose.conformation();
		PDBInfo info = *pose.pdb_info();
		Size ori_nres = conf.size();

		info.attach_to( conf );
		conf.append_polymer_residue_after_seqpos( *ala_rsd, 2, true ); // due to chain ids in input pdb and subsequent generated pose, can only append at position 2
		info.detach_from();

		TS_ASSERT_EQUALS( info.nres(), ori_nres + 1 );
		TS_ASSERT_EQUALS( info.chain( 3 ), PDBInfo::empty_record() );
		TS_ASSERT_EQUALS( info.number( 3 ), 0 );
		TS_ASSERT_EQUALS( info.icode( 3 ), ' ' );
		TS_ASSERT( info.obsolete() );
	}

	/// @brief test Conformation observer behavior on prepend
	void test_PDBInfo_observe_prepend() {
		Conformation & conf = pose.conformation();
		PDBInfo info = *pose.pdb_info();
		Size ori_nres = conf.size();

		info.attach_to( conf );
		conf.prepend_polymer_residue_before_seqpos( *ala_rsd, 3, true ); // due to chain ids in input pdb and subsequent generated pose, can only prepend at position 3
		info.detach_from();

		TS_ASSERT_EQUALS( info.nres(), ori_nres + 1 );
		TS_ASSERT_EQUALS( info.chain( 3 ), PDBInfo::empty_record() );
		TS_ASSERT_EQUALS( info.number( 3 ), 0 );
		TS_ASSERT_EQUALS( info.icode( 3 ), ' ' );
		TS_ASSERT( info.obsolete() );
	}

	/// @brief test Conformation observer behavior on delete
	void test_PDBInfo_observe_delete() {
		Conformation & conf = pose.conformation();
		PDBInfo info = *pose.pdb_info();
		Size ori_nres = conf.size();

		info.attach_to( conf );
		conf.append_polymer_residue_after_seqpos( *ala_rsd, 2, true ); // now 5 residues, need this extra residue because all four residues are jump residues
		conf.delete_polymer_residue( 3 ); // delete the middle residue
		info.detach_from();

		TS_ASSERT_EQUALS( info.nres(), ori_nres );
		TS_ASSERT( info.obsolete() );
	}

	/// @brief test Conformation observer behavior on residue replacement
	void test_PDBInfo_observe_replace() {
		using namespace core::chemical;
		using namespace core::conformation;

		ResidueTypeSetCOP residue_set
			( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
		ResidueOP nt_ala_rsd = ResidueFactory::create_residue( residue_set->name_map( "ALA:NtermProteinFull" ) );

		ConformationOP conf_op = pose.conformation().clone();
		Conformation & conf = *conf_op;
		conf.fold_tree( core::kinematics::FoldTree( 4 ) );
		PDBInfo info = *pose.pdb_info();

		info.attach_to( conf );
		conf.replace_residue( 2, *nt_ala_rsd, true );
		info.detach_from();

		TS_ASSERT_EQUALS( info.natoms( 2 ), nt_ala_rsd->natoms() );
	}

	/// @brief test observer auto-detach when Conformation is destroyed
	void test_PDBInfo_detach_on_destroy() {
		ConformationOP conf( new Conformation( pose.conformation() ) );
		PDBInfo info = *pose.pdb_info();

		info.attach_to( *conf );
		conf = NULL;

		TS_ASSERT( !info.is_observing().lock() );
	}

};
