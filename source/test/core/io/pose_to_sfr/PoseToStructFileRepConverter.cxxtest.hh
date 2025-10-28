// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    test/core/io/pose_to_sfr/PoseToStructFileRepConverter.cxxtest.hh
/// @brief   Test suite for PoseToStructFileRepConverter.
/// @author  Vikram K. Mulligan (vmullig@uw.edu)
/// @details Note that certain functions of PoseToStructFileRepConverter that were originally in
/// file_data are tested in core/io/pose_to_sfr/file_data.cxxtest.hh.

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pdb1ubq.hh>
#include <test/util/pdb1rpb.hh>

// Unit headers
#include <core/io/pose_to_sfr/PoseToStructFileRepConverter.hh>

// Package headers
#include <core/io/StructFileRep.hh>
#include <core/io/AtomInformation.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/membrane/AddMembraneMover.hh>

// Utility headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>

// Basic headers
#include <basic/Tracer.hh>

#include <core/io/HeaderInformation.hh> // AUTO IWYU For HeaderInformation, HeaderInformation::Authors

static basic::Tracer TR("core.io.pose_to_sfr.PoseToStructFileRepConverterTests.cxxtest");

class PoseToStructFileRepConverterTests : public CxxTest::TestSuite {

public:  // Standard methods //////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		core_init_with_additional_options(
			"-write_all_connect_info "
			"-connect_info_cutoff 0.0 "
			"-include_sugars "
			"-output_ligands_as_separate_chains "
			"-INTEGRATION_TEST "
			// Parser doesn't recognize double quotes; using underscore to make my life easier for test.
			"-set_pdb_author J.R.R.TOLKEIN,M.L.KING_JR.,J.W.GIBBS" );
	}

	// Destruction
	void tearDown()
	{}


public:  // Tests /////////////////////////////////////////////////////////////
	// Input Methods //////////////////////////////////////////////////////////

	/// @brief Check that the conect records are properly written out to the StructFileRep.
	void test_conect_records () {
		core::pose::PoseOP testpose( pdb1rpb_poseop() );
		core::io::pose_to_sfr::PoseToStructFileRepConverter pose_to_sfr;
		pose_to_sfr.init_from_pose( *testpose );
		//testpose->dump_pdb("vtemp.pdb");
		// chains is vector0
		for ( core::Size i=0; i<pose_to_sfr.sfr()->chains()[0].size(); ++i ) {
			TS_ASSERT_EQUALS( pose_to_sfr.sfr()->chains()[0][i].serial, i+1); // Atoms should be numbered sequentially.
			// TR << "Atom " << i+1 << ": serial=" << pose_to_sfr.sfr()->chains[1][i-1].serial << " connections=" << pose_to_sfr.sfr()->chains[1][i-1].connected_indices.size() << std::endl;
			TR.flush();
			if ( i == 153 ) {
				TS_ASSERT_EQUALS( pose_to_sfr.sfr()->chains()[0][i].connected_indices.size(), 2 ); //Atom 154 is bound to two atoms.
				TS_ASSERT_EQUALS( pose_to_sfr.sfr()->chains()[0][i].connected_indices[1], 6 ); //Atom 154 is bound to atom 6 (disulfide SG-SG).
				TS_ASSERT_EQUALS( pose_to_sfr.sfr()->chains()[0][i].connected_indices[2], 153 ); //Atom 154 is bound to atom 154 (SG-CB).
			}
		}
	}

	/// @brief Check that the membrane HETATM stuff is properly written out to the StructFileRep.
	void test_membrane_hetatm () {
		// Load in pose from pdb
		core::pose::PoseOP pose( new core::pose::Pose );
		core::import_pose::pose_from_file( *pose, "protocols/membrane/1C3W_TR_A.pdb", core::import_pose::PDB_file );

		//Add membrane information from spanfile:
		protocols::membrane::AddMembraneMoverOP add_memb( new protocols::membrane::AddMembraneMover( "protocols/membrane/1C3W_A.span" ) );
		add_memb->apply( *pose );

		core::io::pose_to_sfr::PoseToStructFileRepConverter pose_to_sfr;
		pose_to_sfr.init_from_pose( *pose );
		pose_to_sfr.grab_membrane_info( *pose, true ); // Re-adds an additional Membrane residue, in addtion to the current in-pose one.

		core::Size const maxchain( pose_to_sfr.sfr()->chains().size() - 1 );
		TS_ASSERT_EQUALS( maxchain, 2 );

		for ( core::Size i=0; i<=maxchain; ++i ) {
			for ( core::Size j=0, jmax=pose_to_sfr.sfr()->chains()[i].size(); j<jmax; ++j ) {
				TR << "chn" << i << " atm" << j << ":\t" << pose_to_sfr.sfr()->chains()[i][j] << std::endl;
			}
		}

		TS_ASSERT_EQUALS(pose_to_sfr.sfr()->chains()[maxchain].size(), 3);
		TS_ASSERT(pose_to_sfr.sfr()->chains()[maxchain][0].isHet);
		TS_ASSERT(pose_to_sfr.sfr()->chains()[maxchain][1].isHet);
		TS_ASSERT(pose_to_sfr.sfr()->chains()[maxchain][2].isHet);
		TS_ASSERT_EQUALS(pose_to_sfr.sfr()->chains()[maxchain][0].serial, 3507);
		TS_ASSERT_EQUALS(pose_to_sfr.sfr()->chains()[maxchain][1].serial, 3508);
		TS_ASSERT_EQUALS(pose_to_sfr.sfr()->chains()[maxchain][2].serial, 3509);
		TS_ASSERT_EQUALS(pose_to_sfr.sfr()->chains()[maxchain][0].name, "THKN");
		TS_ASSERT_EQUALS(pose_to_sfr.sfr()->chains()[maxchain][1].name, "CNTR");
		TS_ASSERT_EQUALS(pose_to_sfr.sfr()->chains()[maxchain][2].name, "NORM");
		TS_ASSERT_EQUALS(pose_to_sfr.sfr()->chains()[maxchain][0].resName, "MEM");
		TS_ASSERT_EQUALS(pose_to_sfr.sfr()->chains()[maxchain][1].resName, "MEM");
		TS_ASSERT_EQUALS(pose_to_sfr.sfr()->chains()[maxchain][2].resName, "MEM");
		TS_ASSERT_EQUALS(pose_to_sfr.sfr()->chains()[maxchain][0].chainID, "C");
		TS_ASSERT_EQUALS(pose_to_sfr.sfr()->chains()[maxchain][1].chainID, "C");
		TS_ASSERT_EQUALS(pose_to_sfr.sfr()->chains()[maxchain][2].chainID, "C");
		TS_ASSERT_EQUALS(pose_to_sfr.sfr()->chains()[maxchain][0].x, 15.3510);
		TS_ASSERT_EQUALS(pose_to_sfr.sfr()->chains()[maxchain][0].y, 0.0 );
		TS_ASSERT_EQUALS(pose_to_sfr.sfr()->chains()[maxchain][0].z, 0.0 );
		TS_ASSERT_EQUALS(pose_to_sfr.sfr()->chains()[maxchain][1].x, 0.0);
		TS_ASSERT_EQUALS(pose_to_sfr.sfr()->chains()[maxchain][1].y, 0.0 );
		TS_ASSERT_EQUALS(pose_to_sfr.sfr()->chains()[maxchain][1].z, 0.0 );
		TS_ASSERT_EQUALS(pose_to_sfr.sfr()->chains()[maxchain][2].x, 0.0);
		TS_ASSERT_EQUALS(pose_to_sfr.sfr()->chains()[maxchain][2].y, 0.0 );
		TS_ASSERT_EQUALS(pose_to_sfr.sfr()->chains()[maxchain][2].z, 15.3510 );
		TS_ASSERT_EQUALS(pose_to_sfr.sfr()->chains()[maxchain][0].element, " H" );
		TS_ASSERT_EQUALS(pose_to_sfr.sfr()->chains()[maxchain][1].element, " H" );
		TS_ASSERT_EQUALS(pose_to_sfr.sfr()->chains()[maxchain][2].element, " H" );
	}

	/// @author  Labonte <JWLabonte@jhu.edu>
	void test_title_section_records () {
		using namespace core;

		pose::PoseOP testpose( pdb1rpb_poseop() );
		io::pose_to_sfr::PoseToStructFileRepConverter pose_to_sfr;
		pose_to_sfr.init_from_pose( *testpose );

		TS_ASSERT_EQUALS( pose_to_sfr.sfr()->header()->authors().size(), 3 );
		TS_ASSERT_EQUALS( pose_to_sfr.sfr()->header()->authors().front(), "J.R.R.TOLKEIN" );
		TS_ASSERT_EQUALS( pose_to_sfr.sfr()->header()->authors().back(), "J.W.GIBBS" );
	}

	/// @brief Check that the foldtree records are properly written out to the StructFileRep.
	void test_foldtree_records () {
		core::pose::PoseOP testpose( pdb1ubq5to13_poseop() );
		core::io::pose_to_sfr::PoseToStructFileRepConverter pose_to_sfr;
		pose_to_sfr.init_from_pose( *testpose );
		pose_to_sfr.grab_foldtree( *testpose, true );

		TS_ASSERT_EQUALS( pose_to_sfr.sfr()->foldtree_string(), "FOLD_TREE  EDGE 1 9 -1 " );
	}

	/// @brief Check that the torsion records are properly written out to the StructFileRep.
	void test_torsion_records () {
		core::pose::PoseOP testpose( new core::pose::Pose );
		core::pose::make_pose_from_sequence( *testpose, "AAAAAA", "fa_standard");

		for ( core::Size i=1, imax=testpose->size(); i<=imax; ++i ) {
			if ( i>1 ) testpose->set_phi(i, -61);
			if ( i<imax ) {
				testpose->set_psi(i, -41);
				testpose->set_omega(i, 180);
			}
		}
		testpose->set_phi(3, -135);
		testpose->set_psi(3, 135);

		core::io::pose_to_sfr::PoseToStructFileRepConverter pose_to_sfr;
		pose_to_sfr.init_from_pose( *testpose );
		pose_to_sfr.grab_torsion_records( *testpose, true );
		if ( TR.visible() ) TR << pose_to_sfr.sfr()->additional_string_output() << std::endl;
		TS_ASSERT_EQUALS( pose_to_sfr.sfr()->additional_string_output(), "REMARK torsions: res    res    chain seq dssp phi psi omega\nREMARK    1    1 1 A L     0.000   -41.000   180.000\nREMARK    2    2 1 A L   -61.000   -41.000   180.000\nREMARK    3    3 1 A L  -135.000   135.000   180.000\nREMARK    4    4 1 A L   -61.000   -41.000   180.000\nREMARK    5    5 1 A L   -61.000   -41.000   180.000\nREMARK    6    6 1 A L   -61.000     0.000     0.000\n" );
	}

	void test_assignment_of_new_chainIDs_for_ligands() {
		using namespace core;
		using namespace io;
		using namespace import_pose;

		TR << "Testing that a .pdb file with a ligand labeled the same as a nearby peptide chain will have the ligands "
			"chain ID reassigned adequately if the option is set to do so." << std::endl;

		pose::Pose pose;
		pose_from_file( pose, "core/io/pose_to_sfr/2_peptides_2_glycan_ligands.pdb", PDB_file);

		// Rosetta reorders the chains by default from how they appear in the .pdb file itself,
		// grouping the three chain As together.
		TS_ASSERT_EQUALS( pose.pdb_info()->chain( 1 ), "A" );  // peptide chain
		TS_ASSERT_EQUALS( pose.pdb_info()->chain( 6 ), "A" );  // oligosaccharide chain
		TS_ASSERT_EQUALS( pose.pdb_info()->chain( 11 ), "A" );  // oligosaccharide chain
		TS_ASSERT_EQUALS( pose.pdb_info()->chain( 16 ), "B" );  // peptide chain

		pose_to_sfr::PoseToStructFileRepConverter pose_to_sfr_converter;
		pose_to_sfr_converter.init_from_pose( pose );
		StructFileRepCOP sfr( pose_to_sfr_converter.sfr() );
		utility::vector0< utility::vector0< AtomInformation > > const & chains( sfr->chains() );

		TS_ASSERT_EQUALS( chains[ 0 ][ 0 ].chainID, "A" );
		TS_ASSERT_EQUALS( chains[ 1 ][ 0 ].chainID, "C" );
		TS_ASSERT_EQUALS( chains[ 2 ][ 0 ].chainID, "D" );
		TS_ASSERT_EQUALS( chains[ 3 ][ 0 ].chainID, "B" );
	}

	void test_ssbond_numbering() {
		using namespace core;
		using namespace io;

		core::pose::PoseOP testpose( pdb1rpb_poseop() );
		// Make sure we shift things.
		testpose->delete_polymer_residue( 5 );
		testpose->delete_polymer_residue( 4 );
		testpose->delete_polymer_residue( 3 );

		StructFileRepOptions opts;
		opts.set_renumber_pdb( true );
		core::io::pose_to_sfr::PoseToStructFileRepConverter pose_to_sfr(opts);
		pose_to_sfr.init_from_pose( *testpose );
		StructFileRepOP sfr = pose_to_sfr.sfr();

		// In original, CYS 1A is bound to CYS 13A & CYS 7A to CYS 19A
		// In the delete+renumber, it should be 1A-10A && 4A-16A
		TS_ASSERT_EQUALS( sfr->ssbond_map().size(), 2 );
		for ( auto const & entry: sfr->ssbond_map() ) {
			TS_ASSERT_EQUALS( entry.second.size(), 1 );
			SSBondInformation const & ssbi = entry.second[1];
			if ( ssbi.resSeq1 == 1 ) {
				TS_ASSERT_EQUALS( ssbi.resSeq1, 1 );
				TS_ASSERT_EQUALS( ssbi.chainID1, "A" );
				TS_ASSERT_EQUALS( ssbi.resSeq2, 10 );
				TS_ASSERT_EQUALS( ssbi.chainID2, "A" );
			} else {
				TS_ASSERT_EQUALS( ssbi.resSeq1, 4 );
				TS_ASSERT_EQUALS( ssbi.chainID1, "A" );
				TS_ASSERT_EQUALS( ssbi.resSeq2, 16 );
				TS_ASSERT_EQUALS( ssbi.chainID2, "A" );
			}
		}

		// Double check that the Atom information has the appropriate settings
		for ( auto const & chain: sfr->chains() ) {
			for ( auto const & ai: chain ) {
				if ( ai.resSeq == 1 || ai.resSeq == 10 || ai.resSeq == 4 || ai.resSeq == 16 ) {
					TS_ASSERT_EQUALS( ai.resName, "CYS" );
				}
			}
		}
	}

	//SEE ALSO StructFileRep.cxxtest.hh; unclear if test_generate_secondary_structure_informations belongs here or there

};  // class StructFileRepTests
