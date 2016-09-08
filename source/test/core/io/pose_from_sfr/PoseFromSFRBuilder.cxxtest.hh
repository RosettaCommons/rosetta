// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/pose_builder/PoseFromSFRBuilder.cxxtest.hh
/// @brief  test suite for classes associated with core::io::pose_builder::PoseFromSFRBuilder
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Unit headers
#include <core/io/pose_from_sfr/PoseFromSFRBuilder.hh>
#include <core/io/HeaderInformation.hh>

// Program headers
#include <core/io/StructFileRep.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/io/pdb/pdb_reader.hh>
#include <core/io/StructFileReaderOptions.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

// Basic headers
#include <basic/Tracer.hh>

static basic::Tracer TR("core.io.pose_from_sfr.PoseFromSFRBuilder.cxxtest");

using namespace core;
using namespace core::io;
using namespace ObjexxFCL;
using namespace core::io::pose_from_sfr;

class PoseFromSFRBuilderTests : public CxxTest::TestSuite
{

public:
	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	void test_build_trpcage_pdb() {
		pose::Pose pose;
		core::io::StructFileReaderOptions options;
		core::io::StructFileRep sfr = core::io::pdb::create_sfr_from_pdb_file_contents( trp_cage_ideal(), options );
		chemical::ResidueTypeSetCOP residue_set
			( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ) );
		PoseFromSFRBuilder pb( residue_set, options );
		pb.build_pose( sfr, pose );

		TS_ASSERT( pose.size() == 20 );
		if ( pose.size() != 20 ) return;

		TS_ASSERT( pose.residue_type(  1 ).has_variant_type( chemical::LOWER_TERMINUS_VARIANT ));
		for ( Size ii = 2; ii <= 19; ++ii ) {
			TS_ASSERT( ! pose.residue_type( ii ).has_variant_type( chemical::LOWER_TERMINUS_VARIANT ));
			TS_ASSERT( ! pose.residue_type( ii ).has_variant_type( chemical::UPPER_TERMINUS_VARIANT ));
		}
		TS_ASSERT( pose.residue_type( 20 ).has_variant_type( chemical::UPPER_TERMINUS_VARIANT ));
		//pose.dump_pdb( "test_pose_builder.pdb" );
	}

	void test_read_in_acylated_nterm_residue() {
		core::io::StructFileReaderOptions options;
		std::string ace_as_own_residue_pdb = utility::file_contents( "core/io/pose_from_sfr/ace_from_1VDN.pdb" );
		core::io::StructFileRep sfr = core::io::pdb::create_sfr_from_pdb_file_contents( ace_as_own_residue_pdb, options );

		chemical::ResidueTypeSetCOP residue_set
			( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ) );
		PoseFromSFRBuilder builder( residue_set, options );
		pose::Pose pose;
		builder.build_pose( sfr, pose );

		TS_ASSERT_EQUALS( pose.size(), 4 );
		TS_ASSERT_EQUALS( pose.residue(1).aa(), chemical::aa_ala );
		TS_ASSERT_EQUALS( pose.residue(2).aa(), chemical::aa_ala );
		TS_ASSERT_EQUALS( pose.residue(3).aa(), chemical::aa_pro );
		TS_ASSERT_EQUALS( pose.residue(4).aa(), chemical::aa_ala );
		TS_ASSERT( pose.residue_type(1).has_property( "ACETYLATED_NTERMINUS" ));
		if ( ! pose.residue_type(1).has_property( "ACETYLATED_NTERMINUS" ) ) return;

		id::AtomID_Mask const & missing = builder.missing_atoms();
		TS_ASSERT( ! missing[ id::AtomID( pose.residue_type(1).atom_index( "CO" ),  1 ) ] );
		TS_ASSERT( ! missing[ id::AtomID( pose.residue_type(1).atom_index( "OP1" ), 1 ) ] );
		TS_ASSERT(   missing[ id::AtomID( pose.residue_type(1).atom_index( "CP2" ), 1 ) ] );

		//pose.dump_file( "fill_ace_missing_atoms.pdb" );
	}

	void test_read_in_res_both_acetylated_and_cterm() {
		// This isn't really a test of the PoseFromSFRBuilder so much as it is a test of the
		// AcetylatedProteinNTerm and MethylatedProteinCTerm patches and of the ResidueTypeFinder

		core::io::StructFileReaderOptions options;
		std::string ace_and_meth_ala_pdb = utility::file_contents( "core/io/pose_from_sfr/acetylated_and_methylated_alanine.pdb" );
		core::io::StructFileRep sfr = core::io::pdb::create_sfr_from_pdb_file_contents( ace_and_meth_ala_pdb, options );

		chemical::ResidueTypeSetCOP residue_set
			( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ) );
		PoseFromSFRBuilder builder( residue_set, options );
		pose::Pose pose;
		builder.build_pose( sfr, pose );

		TS_ASSERT_EQUALS( pose.size(), 1 );
		TS_ASSERT_EQUALS( pose.residue(1).aa(), chemical::aa_ala );
		TS_ASSERT( pose.residue_type(1).has_property( "ACETYLATED_NTERMINUS" ));
		TS_ASSERT( pose.residue_type(1).has_property( "METHYLATED_CTERMINUS" ));
	}

};
