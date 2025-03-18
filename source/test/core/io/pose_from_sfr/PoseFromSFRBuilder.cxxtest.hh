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
#include <test/util/pdb1lnt.hh>

// Unit headers
#include <core/io/pose_from_sfr/PoseFromSFRBuilder.hh>

// Program headers
#include <core/io/StructFileRep.hh>
#include <core/io/pdb/pdb_reader.hh>
#include <core/io/StructFileReaderOptions.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/import_pose_options.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>

static basic::Tracer TR( "core.io.pose_from_sfr.PoseFromSFRBuilder.cxxtest" );

using namespace core;
using namespace core::io;
using namespace ObjexxFCL;
using namespace core::io::pose_from_sfr;

class PoseFromSFRBuilderTests : public CxxTest::TestSuite
{

public:
	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-include_sugars" );
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

		TR << "Pose residue 1: " << pose.residue_type(1).name() << std::endl;

		TS_ASSERT_EQUALS( pose.size(), 1 );
		TS_ASSERT_EQUALS( pose.residue(1).aa(), chemical::aa_ala );
		TS_ASSERT( pose.residue_type(1).has_property( "ACETYLATED_NTERMINUS" ));
		TS_ASSERT( pose.residue_type(1).has_property( "METHYLATED_CTERMINUS" ));
	}

	// Confirm that pass_1_merge_and_split_residues_as_necessary() correctly splits a splitable residue.
	/// @author  Labonte <JWLabonte@jhu.edu>
	void test_splitting_of_residues() {
		using namespace std;
		using namespace chemical;

		TR << "Testing that the PoseFromSFRBuilder can split residues." << endl;

		StructFileReaderOptions options;
		string const pdb_contents(  // Lines from PDB #5BJZ
			"ATOM     36  N   GLY A   5      -9.957 -28.821 189.708  1.00 36.49           N  \n"
			"ATOM     37  CA  GLY A   5     -10.536 -29.933 188.980  1.00 33.76           C  \n"
			"ATOM     38  C   GLY A   5     -11.743 -29.607 188.130  1.00 32.07           C  \n"
			"ATOM     39  O   GLY A   5     -12.472 -30.526 187.740  1.00 30.93           O  \n"
			"HETATM12532  C1  MAL B 401     -16.780  -1.670 175.019  1.00 22.31           C  \n"
			"HETATM12533  C2  MAL B 401     -15.458  -1.864 175.750  1.00 20.10           C  \n"
			"HETATM12534  C3  MAL B 401     -15.725  -1.759 177.260  1.00 22.44           C  \n"
			"HETATM12535  C4  MAL B 401     -16.356  -0.412 177.578  1.00 23.54           C  \n"
			"HETATM12536  C5  MAL B 401     -17.612  -0.218 176.728  1.00 21.74           C  \n"
			"HETATM12537  C6  MAL B 401     -18.181   1.184 176.949  1.00 24.56           C  \n"
			"HETATM12538  O1  MAL B 401     -17.679  -2.694 175.383  1.00 23.77           O  \n"
			"HETATM12539  O2  MAL B 401     -14.947  -3.131 175.441  1.00 19.29           O  \n"
			"HETATM12540  O3  MAL B 401     -14.531  -1.936 177.997  1.00 20.90           O  \n"
			"HETATM12541  O4  MAL B 401     -16.698  -0.306 178.948  1.00 22.31           O  \n"
			"HETATM12542  O5  MAL B 401     -17.323  -0.414 175.357  1.00 21.30           O  \n"
			"HETATM12543  O6  MAL B 401     -19.358   1.387 176.189  1.00 24.58           O  \n"
			"HETATM12544  C1' MAL B 401     -20.164  -5.437 173.448  1.00 27.37           C  \n"
			"HETATM12545  C2' MAL B 401     -18.684  -5.407 173.044  1.00 26.37           C  \n"
			"HETATM12546  C3' MAL B 401     -17.785  -4.643 174.017  1.00 24.36           C  \n"
			"HETATM12547  C4' MAL B 401     -18.418  -3.307 174.341  1.00 26.47           C  \n"
			"HETATM12548  C5' MAL B 401     -19.887  -3.472 174.771  1.00 25.02           C  \n"
			"HETATM12549  C6' MAL B 401     -20.564  -2.126 174.943  1.00 21.11           C  \n"
			"HETATM12550  O1' MAL B 401     -20.398  -6.331 174.513  1.00 29.23           O  \n"
			"HETATM12551  O2' MAL B 401     -18.240  -6.733 172.968  1.00 32.21           O  \n"
			"HETATM12552  O3' MAL B 401     -16.481  -4.449 173.475  1.00 26.14           O  \n"
			"HETATM12553  O5' MAL B 401     -20.646  -4.162 173.796  1.00 24.94           O  \n"
			"HETATM12554  O6' MAL B 401     -21.712  -2.306 175.755  1.00 27.37           O  \n"
			"ATOM   2911  N   GLY C  11     -54.011   1.690 100.082  1.00 29.22           N  \n"
			"ATOM   2912  CA  GLY C  11     -54.910   0.571  99.922  1.00 26.67           C  \n"
			"ATOM   2913  C   GLY C  11     -54.864  -0.464 101.023  1.00 31.14           C  \n"
			"ATOM   2914  O   GLY C  11     -55.627  -1.435 100.960  1.00 35.10           O  \n" );
		StructFileRep const sfr( pdb::create_sfr_from_pdb_file_contents( pdb_contents, options ) );

		ResidueTypeSetCOP residue_set( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
		PoseFromSFRBuilder builder( residue_set, options );
		pose::Pose pose;
		builder.build_pose( sfr, pose );

		TS_ASSERT_EQUALS( pose.size(), 4 );

		pose::PDBInfoCOP info( pose.pdb_info() );

		TS_ASSERT_EQUALS( pose.residue( 1 ).name3(), "GLY" );
		TS_ASSERT_EQUALS( info->chain( 1 ), 'A' );
		TS_ASSERT_EQUALS( info->number( 1 ), 5 );
		TS_ASSERT_EQUALS( pose.residue( 2 ).name3(), "Glc" );
		TS_ASSERT_EQUALS( info->chain( 2 ), 'B' );
		TS_ASSERT_EQUALS( info->number( 2 ), 401 );
		TS_ASSERT_EQUALS( info->icode( 2 ), 'A' );
		TS_ASSERT_EQUALS( pose.residue( 3 ).name3(), "Glc" );
		TS_ASSERT_EQUALS( info->chain( 3 ), 'B' );
		TS_ASSERT_EQUALS( info->number( 3 ), 401 );
		TS_ASSERT_EQUALS( info->icode( 3 ), 'B' );
		TS_ASSERT_EQUALS( pose.residue( 4 ).name3(), "GLY" );
		TS_ASSERT_EQUALS( info->chain( 4 ), 'C' );
		TS_ASSERT_EQUALS( info->number( 4 ), 11 );
	}

	void test_1rgr() {
		core::io::StructFileReaderOptions options;
		std::string ace_as_own_residue_pdb = utility::file_contents( "core/io/pose_from_sfr/1rgr_apo.pdb" );
		core::io::StructFileRep sfr = core::io::pdb::create_sfr_from_pdb_file_contents( ace_as_own_residue_pdb, options );

		chemical::ResidueTypeSetCOP residue_set
			( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ) );
		PoseFromSFRBuilder builder( residue_set, options );
		pose::Pose pose;
		builder.build_pose( sfr, pose );

		// This pose has the following chemical connectivity:
		// Residues 1-6 are connected polymerically as normal.
		// Residue 3's sidechain is connected to residue 7's Cterm
		// Residue 5's sidechain is connected to residue 7's Nterm
		TS_ASSERT_EQUALS( pose.residue(1).n_current_residue_connections(), 1 );
		TS_ASSERT_EQUALS( pose.residue(2).n_current_residue_connections(), 2 );
		TS_ASSERT_EQUALS( pose.residue(3).n_current_residue_connections(), 3 );
		TS_ASSERT_EQUALS( pose.residue(4).n_current_residue_connections(), 2 );
		TS_ASSERT_EQUALS( pose.residue(5).n_current_residue_connections(), 3 );
		TS_ASSERT_EQUALS( pose.residue(6).n_current_residue_connections(), 1 );
		TS_ASSERT_EQUALS( pose.residue(7).n_current_residue_connections(), 2 );


		// But enough with chemical connectivity. What about a reasonable foldtree?
		// When does that happen?
		// OK, this isn't properly a test of the PoseFromSFRBuilder anymore.
		core::pose::set_reasonable_fold_tree( pose );
		TS_ASSERT_EQUALS( pose.fold_tree().num_jump(), 0 );

	}

	/// @brief This tests the ability for Rosetta to read in a cyclic peptide that has a LINK
	/// record connecting the termini.  It should automatically set up a covalent bond and
	/// cyclization variants.
	void test_read_cyclic_peptide_with_link_record() {
		core::pose::PoseOP pose( core::import_pose::pose_from_file( "core/io/pose_from_sfr/cyclic_pep_with_link.pdb", false ) );
		TS_ASSERT( pose->residue_type(1).has_variant_type( core::chemical::CUTPOINT_UPPER ) );
		TS_ASSERT( pose->residue_type(pose->total_residue()).has_variant_type( core::chemical::CUTPOINT_LOWER ) );
		TS_ASSERT( pose->residue(1).connected_residue_at_lower() == pose->total_residue() );
		TS_ASSERT( pose->residue(pose->total_residue()).connected_residue_at_upper() == 1 );
	}

	/// @brief This tests for proper chain separation when you have a cyclic peptide after a regular one.
	/// There was a bug where the peptide cycle was counted as "part of" the previous chain,
	/// due to the lack of being "terminal"
	void test_read_cyclic_peptide_chain_ending() {
		core::pose::PoseOP pose( core::import_pose::pose_from_file( "core/io/pose_from_sfr/cycle_chain_test.pdb", false ) );
		TS_ASSERT_EQUALS( pose->size(), 5 );
		TR << "cycle_chain_test.pdb sequence: " << pose->annotated_sequence() << std::endl;
		TS_ASSERT( pose->residue(2).is_upper_terminus() );
		TS_ASSERT( !pose->residue(3).is_lower_terminus() );
		TS_ASSERT( !pose->residue(5).is_upper_terminus() );
		TS_ASSERT_EQUALS( pose->residue(3).connected_residue_at_lower(), 5 );
		TS_ASSERT_EQUALS( pose->residue(5).connected_residue_at_upper(), 3 );
		TS_ASSERT_EQUALS( pose->num_chains(), 2 );
	}

	/// @brief Test that we can correctly read poses without adding termini to user-specified chains.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	void test_read_pose_without_terminal_variant_types() {
		core::import_pose::ImportPoseOptions options1, options2, options3, options4;
		options1.set_pack_missing_sidechains(false);
		options2.set_pack_missing_sidechains(false);
		options3.set_pack_missing_sidechains(false);
		options4.set_pack_missing_sidechains(false);
		options2.set_check_if_residues_are_Ntermini( "" );
		options2.set_check_if_residues_are_Ctermini( "" );
		options3.set_check_if_residues_are_Ntermini( "B" );
		options3.set_check_if_residues_are_Ctermini( "B" );
		options4.set_check_if_residues_are_Ntermini( "AC" );
		options4.set_check_if_residues_are_Ctermini( "AC" );

		core::pose::PoseOP pose1( core::import_pose::pose_from_file( "core/io/pose_from_sfr/1IJ3_cleaned.pdb" , options1, false, core::import_pose::PDB_file  ) );
		core::pose::PoseOP pose2( core::import_pose::pose_from_file( "core/io/pose_from_sfr/1IJ3_cleaned.pdb" , options2, false, core::import_pose::PDB_file  ) );
		core::pose::PoseOP pose3( core::import_pose::pose_from_file( "core/io/pose_from_sfr/1IJ3_cleaned.pdb" , options3, false, core::import_pose::PDB_file  ) );
		core::pose::PoseOP pose4( core::import_pose::pose_from_file( "core/io/pose_from_sfr/1IJ3_cleaned.pdb" , options4, false, core::import_pose::PDB_file  ) );

		core::Size const res1A(1), resNA(31), res1B(32), resNB(62), res1C(63), resNC(93); //Start and end indices of all chains in the test structure.

		//Pose 1 should have termini.
		TS_ASSERT( pose1->residue_type(res1A).is_lower_terminus() );
		TS_ASSERT( pose1->residue_type(resNA).is_upper_terminus() );
		TS_ASSERT( pose1->residue_type(res1B).is_lower_terminus() );
		TS_ASSERT( pose1->residue_type(resNB).is_upper_terminus() );
		TS_ASSERT( pose1->residue_type(res1C).is_lower_terminus() );
		TS_ASSERT( pose1->residue_type(resNC).is_upper_terminus() );

		//Pose 2 should have no termini.
		TS_ASSERT( !pose2->residue_type(res1A).is_lower_terminus() );
		TS_ASSERT( !pose2->residue_type(resNA).is_upper_terminus() );
		TS_ASSERT( !pose2->residue_type(res1B).is_lower_terminus() );
		TS_ASSERT( !pose2->residue_type(resNB).is_upper_terminus() );
		TS_ASSERT( !pose2->residue_type(res1C).is_lower_terminus() );
		TS_ASSERT( !pose2->residue_type(resNC).is_upper_terminus() );

		//Pose 3 should have termini ONLY on chain B.
		TS_ASSERT( !pose3->residue_type(res1A).is_lower_terminus() );
		TS_ASSERT( !pose3->residue_type(resNA).is_upper_terminus() );
		TS_ASSERT( pose3->residue_type(res1B).is_lower_terminus() );
		TS_ASSERT( pose3->residue_type(resNB).is_upper_terminus() );
		TS_ASSERT( !pose3->residue_type(res1C).is_lower_terminus() );
		TS_ASSERT( !pose3->residue_type(resNC).is_upper_terminus() );

		//Pose 4 should have termini ONLY on chains A and C.
		TS_ASSERT( pose4->residue_type(res1A).is_lower_terminus() );
		TS_ASSERT( pose4->residue_type(resNA).is_upper_terminus() );
		TS_ASSERT( !pose4->residue_type(res1B).is_lower_terminus() );
		TS_ASSERT( !pose4->residue_type(resNB).is_upper_terminus() );
		TS_ASSERT( pose4->residue_type(res1C).is_lower_terminus() );
		TS_ASSERT( pose4->residue_type(resNC).is_upper_terminus() );
	}

	/// @brief utility function for pose comparison
	void compare_poses( core::pose::Pose const & plain, core::pose::Pose const & fast, std::string const & filename ) {
		// Basic check to see if we got poses of the same size
		TSM_ASSERT_EQUALS( filename, plain.size(), fast.size() );
		if ( plain.size() != fast.size() ) {
			TR.Error << filename << " PLAIN " << plain.annotated_sequence() << std::endl;
			TR.Error << filename << " FAST " << fast.annotated_sequence() << std::endl;
			return;
		}

		// Check residue typing
		bool failed = false;
		for ( core::Size ii(1); ii <= fast.size(); ++ii ) {
			if ( plain.residue_type(ii).name() != fast.residue_type(ii).name() ) {
				TR.Error << filename << " Residue " << ii << std::endl;
				TSM_ASSERT_EQUALS( filename, plain.residue_type(ii).name(), fast.residue_type(ii).name() );
				failed = true;
			}
		}
		if ( failed ) { return; }

		// Commented out for speed purposes
		//   // Double check atom assignments
		//   for( core::Size ii(1); ii <= fast.size(); ++ii ) {
		//    core::conformation::Residue const & pres = plain.residue(ii);
		//    core::conformation::Residue const & fres = fast.residue(ii);
		//    for ( core::Size aa(1); aa <= fres.natoms(); ++aa ) {
		//     if ( pres.is_virtual( aa ) ) { continue; } // Not necessarily the same
		//     if (  pres.xyz(aa) != fres.xyz(aa) ) {
		//      TR.Error << filename << " Residue " << ii << " " << pres.type().name() << " atom " << pres.atom_name(aa) << std::endl;
		//      TSM_ASSERT_EQUALS( filename , pres.xyz(aa), fres.xyz(aa) );
		//      failed = true;
		//     }
		//    }
		//   }
		//   if ( failed ) { return; }
	}

	void load_fast_and_slow( std::string const & file_contents, std::string const & filename, chemical::ResidueTypeSetCOP residue_set ) {
		StructFileReaderOptions plain_opt;
		plain_opt.set_fast_restyping(false);
		StructFileReaderOptions fast_opt;
		fast_opt.set_fast_restyping(true);

		TR << "==================== Testing " << filename << " ========================" << std::endl;
		core::pose::Pose plain;
		core::pose::Pose fast;

		core::io::StructFileRep sfr_plain = core::io::pdb::create_sfr_from_pdb_file_contents( file_contents, plain_opt );
		sfr_plain.filename() = filename;
		PoseFromSFRBuilder pb_plain( residue_set, plain_opt );
		pb_plain.build_pose( sfr_plain, plain );

		TR << filename << " PLAIN " << plain.annotated_sequence() << std::endl;
		TR << "-------- Loading fast ------------- " << std::endl;

		core::io::StructFileRep sfr_fast = core::io::pdb::create_sfr_from_pdb_file_contents( file_contents, fast_opt );
		sfr_fast.filename() = filename;
		PoseFromSFRBuilder pb_fast( residue_set, fast_opt );
		pb_fast.build_pose( sfr_fast, fast );

		compare_poses( plain, fast, filename );
	}

	void test_quick_and_dirty_restyping() {
		using namespace core::import_pose;
		using namespace core::io;

		utility::vector1< std::string > fa_files_to_test = {
			"core/io/1QYS.pdb", // Standard PDB loading.
			"core/io/test_in.pdb", // Issue with HIS/HIS_D calling.
			"protocols/sparta/2kywA.pdb", // HIS/HIS_D calling with both protons
			"core/io/daa.pdb", // D-AA
			"core/conformation/4gatA.pdb", // Metal ions & DNA
			"protocols/stepwise/modeler/align/scaff_subset_stepwise_input_seq_1.pdb", // RNA
			"core/conformation/symmetry/mtest1.pdb", // VRT and INV_VRT
			};
		//// Explicitly not supported for the Q&D approach
		//   "core/io/mmtf/crosslinkermover_octahedral_s2_symm_S_0008.pdb", // VirtualMetalConjugation patch due to presence of VM1 atom.
		//   "core/io/two_lipids.pdb", // Lipid names aren't their three letter codes.
		//      "protocols/match/E1cb_carbaryl_1his_oxy_1bb_10_2.pdb", // Needs params, and params don't match with three letter code
		//   "protocols/simple_moves/oop/oop_test.pdb", // OOP -- needs atom annotations
		//   "core/chemical/carbohydrates/gp120_2glycans_man5.pdb", // Glycans w/ HETNAM without patching info
		//   "core/chemical/carbohydrates/GlcN.pdb" // No patching info
		//   "core/chemical/carbohydrates/alpha-L-Fucp-_1-6_-D-GlcpNAc-_1-4_-D-GlcpNAc.pdb", // No patching info
		//   "core/io/pose_from_sfr/acetylated_and_methylated_alanine.pdb", // Terminus alterations -- need patching information
		//"core/io/5FYL.pdb", // Glycans -- Strange issue with ignoring BMA in plain

		utility::vector1< std::string > cen_files_to_test = {
			"devel/znhash/1EER_cn4_0043.pdb",  // Centroid, VRT
			};

		chemical::ResidueTypeSetCOP fa_residue_set( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ) );

		for ( std::string const & filename: fa_files_to_test ) {
			std::string const & file_contents = utility::file_contents( filename );

			load_fast_and_slow( file_contents, filename, fa_residue_set );
		}

		// load_fast_and_slow( pdb1lnt(), "pdb1lnt", fa_residue_set ); //Residue 1 -- "URA:LowerRNA" != "URA:LowerRNA:Virtual_Phosphate"

		chemical::ResidueTypeSetCOP cen_residue_set( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID ) );
		for ( std::string const & filename: cen_files_to_test ) {
			std::string const & file_contents = utility::file_contents( filename );

			load_fast_and_slow( file_contents, filename, cen_residue_set );
		}
	}

};
