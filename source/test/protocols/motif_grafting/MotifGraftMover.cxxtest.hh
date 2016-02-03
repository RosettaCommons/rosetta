// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/motif_grafting/MotifGraftMover.cxxtest.hh
/// @brief Test suite for protocols/motif_grafting/moves/MotifGraftMover
//
// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <util/pose_funcs.hh>

// Project headers
#include <basic/Tracer.hh>
#include <basic/options/keys/sicdock.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>

//Include Rosetta Core Stuff
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/VariantType.hh>

//Include Rosetta numeric
#include <numeric/xyz.functions.hh>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

//My own class headers
#include <protocols/motif_grafting/movers/MotifGraftMover.hh>


static basic::Tracer TR("protocols.motif_grafting.MotifGraftMover.cxxtest.hh");

namespace {
class MotifGraftMoverTests : public CxxTest::TestSuite {

public:

	// Shared set up goes here
	void setUp()
	{
		core_init_with_additional_options("");
	}

	// Shared finalization goes here.
	void tearDown() {

	}

	void test_self_graft_recovery_oneFragment_using_fixBBmatching()
	{
		using namespace protocols::motif_grafting::movers;

		TR << "Start Test 1" << std::endl;

		//Create the mover
		TR << "Creating mover" << std::endl;
		MotifGraftMover mgm;
		TR << "Intializing mover" << std::endl;

		// initialize the parameters
		std::string s_context     = "protocols/motif_grafting/2nm1_context.pdb";       //context structure
		std::string s_motif      = "protocols/motif_grafting/1zx3a_as_motif_1fragment.pdb"; //motif with one fragment
		core::Real r_RMSD      = 0.5;               //RMSD tol
		core::Real r_NC_RMSD     = 0.5;               //NC points RMSD tol
		core::Size  s_ClashTol    = 5;                //Clash Score tol
		std::string s_comb_d     = "0:0";              //Combinatory delta
		std::string s_max_frag_d   = "0:0";              //Max Fragment relpacement delta
		std::string s_clash_res    = "GLY";              //Clash score test residue
		std::string s_hotspots    = "";                //hotspots
		bool  b_full_bb_alignment  = false;              //full_bb_alignment
		bool  b_allow_ind_frag_alig = false;              //allow_independent_alignment_per_fragment
		bool  b_only_copy_sidechain = false;              //only copy sidechains for hotspots (opposed to graft)
		bool  b_only_if_Npoint_eq  = false;              //graft only if Npoint aminoacid identities are equal
		bool  b_only_if_Cpoint_eq  = false;              //graft only if Cpoint aminoacid identities are equal
		bool  b_revert_to_native  = false;              //revert to the native scafold sequence after grafting
		bool  b_allow_repeat    = false;              //Allow repeat same graft

		mgm.init_parameters(
			s_context,
			s_motif,
			r_RMSD,
			r_NC_RMSD,
			s_ClashTol,
			s_comb_d,
			s_max_frag_d,
			s_clash_res,
			s_hotspots,
			b_full_bb_alignment,
			b_allow_ind_frag_alig,
			b_only_copy_sidechain,
			b_only_if_Npoint_eq,
			b_only_if_Cpoint_eq,
			b_revert_to_native,
			b_allow_repeat);

		//Import the test scaffold
		TR << "Importing Pose" << std::endl;
		core::pose::Pose testPose;
		core::import_pose::pose_from_file( testPose, "protocols/motif_grafting/1zx3a_translated_and_rotated.pdb" , core::import_pose::PDB_file);

		//Get the original sequence
		std::string seq_motif_before=testPose.chain_sequence(1);

		TR << "Applying mover" << std::endl;
		//Apply the mover
		mgm.apply(testPose);

		//get the sequence of the epigraft
		std::string seq_motif_after=testPose.chain_sequence(2);

		//Check that the sequences are equal:
		TS_ASSERT(seq_motif_before.size()==seq_motif_after.size());
		for ( core::Size i=0; i < seq_motif_before.size(); ++i ) {
			//TR << seq_motif_before.at(i) << " : " << seq_motif_after.at(i) << std::endl;
			TS_ASSERT( seq_motif_before.at(i) == seq_motif_after.at(i) );
		}

		/////Import the target scaffold "Real Answer"
		core::pose::Pose challengePose;
		core::import_pose::pose_from_file( challengePose, "protocols/motif_grafting/2nm1_1zx3a_hybrid.pdb" , core::import_pose::PDB_file);


		/////Check that "Answer" and the "Real Answer" sequences are equal size:
		TS_ASSERT( testPose.total_residue() == challengePose.total_residue() );

		/////Check that the atomic positions of the "Answer" and the "Real Answer" are equal
		for ( core::Size i = 1; i <= testPose.total_residue(); ++i ) {
			//core::conformation::ResidueOP resA = new core::conformation::Residue (testPose.residue( i ) );
			//core::conformation::ResidueOP resB = new core::conformation::Residue (challengePose.residue( i ) );
			for ( core::Size j = 1; j <= testPose.residue_type(i).natoms(); ++j ) {
				//This is because of the hack on the 1st hydrogen (1H) of the first motif
				if ( !testPose.residue(i).atom_is_hydrogen(j) ) {
					core::id::AtomID id( j, i );
					//TR << testPose.xyz(id)(1) << " : " << challengePose.xyz(id)(1) << std::endl;
					TS_ASSERT_DELTA( testPose.xyz(id), challengePose.xyz(id), 0.001 );
				}
			}
		}
		TR << "End Test" << std::endl;
	}


	void test_self_graft_recovery_twoFragment_using_fixBBmatching()
	{
		using namespace protocols::motif_grafting::movers;

		TR << "Start Test 2" << std::endl;

		//Create the mover
		TR << "Creating mover" << std::endl;
		MotifGraftMover mgm;
		TR << "Intializing mover" << std::endl;

		// initialize the parameters
		std::string s_context     = "protocols/motif_grafting/2nm1_context.pdb";       //context structure
		std::string s_motif      = "protocols/motif_grafting/1zx3a_as_motif_2fragment.pdb"; //motif with one fragment
		core::Real r_RMSD      = 0.5;               //RMSD tol
		core::Real r_NC_RMSD     = 0.5;               //NC points RMSD tol
		core::Size  s_ClashTol    = 5;                //Clash Score tol
		std::string s_comb_d     = "0:0,0:0";            //Combinatory delta
		std::string s_max_frag_d   = "0:0,0:0";            //Max Fragment relpacement delta
		std::string s_clash_res    = "GLY";              //Clash score test residue
		std::string s_hotspots    = "";                //hotspots
		bool   b_full_bb_alignment  = true;               //full_bb_alignment
		bool   b_allow_ind_frag_alig = false;              //allow_independent_alignment_per_fragment
		bool   b_only_copy_sidechain = false;              //only copy sidechains for hotspots (opposed to graft)
		bool  b_only_if_Npoint_eq  = false;              //graft only if Npoint aminoacid identities are equal
		bool  b_only_if_Cpoint_eq  = false;              //graft only if Cpoint aminoacid identities are equal
		bool   b_revert_to_native  = true;               //revert to the native scafold sequence after grafting
		bool   b_allow_repeat    = false;              //Allow repeat same graft

		mgm.init_parameters(
			s_context,
			s_motif,
			r_RMSD,
			r_NC_RMSD,
			s_ClashTol,
			s_comb_d,
			s_max_frag_d,
			s_clash_res,
			s_hotspots,
			b_full_bb_alignment,
			b_allow_ind_frag_alig,
			b_only_copy_sidechain,
			b_only_if_Npoint_eq,
			b_only_if_Cpoint_eq,
			b_revert_to_native,
			b_allow_repeat);

		//Import the test scaffold
		TR << "Importing Pose" << std::endl;
		core::pose::Pose testPose;
		core::import_pose::pose_from_file( testPose, "protocols/motif_grafting/1zx3a_translated_and_rotated.pdb" , core::import_pose::PDB_file);

		//Get the original sequence
		std::string seq_motif_before=testPose.chain_sequence(1);

		TR << "Applying mover" << std::endl;
		//Apply the mover
		mgm.apply(testPose);

		TR << "Mover Finished" << std::endl;

		//get the sequence of the epigraft
		std::string seq_motif_after=testPose.chain_sequence(2);

		//Check that the sequences are equal:
		TS_ASSERT(seq_motif_before.size()==seq_motif_after.size());
		for ( core::Size i=0; i < seq_motif_before.size(); ++i ) {
			//TR << seq_motif_before.at(i) << " : " << seq_motif_after.at(i) << std::endl;
			TS_ASSERT( seq_motif_before.at(i) == seq_motif_after.at(i) );
		}

		/////Import the target scaffold "Real Answer"
		core::pose::Pose challengePose;
		core::import_pose::pose_from_file( challengePose, "protocols/motif_grafting/2nm1_1zx3a_hybrid.pdb" , core::import_pose::PDB_file);


		/////Check that "Answer" and the "Real Answer" sequences are equal size:
		TS_ASSERT( testPose.total_residue() == challengePose.total_residue() );

		/////Check that the atomic positions of the "Answer" and the "Real Answer" are equal
		for ( core::Size i = 1; i <= testPose.total_residue(); ++i ) {
			//core::conformation::ResidueOP resA = new core::conformation::Residue (testPose.residue( i ) );
			//core::conformation::ResidueOP resB = new core::conformation::Residue (challengePose.residue( i ) );
			for ( core::Size j = 1; j <= testPose.residue_type(i).natoms(); ++j ) {
				//This is because of the hack on the 1st hydrogen (1H) of the first motif
				if ( !testPose.residue(i).atom_is_hydrogen(j) ) {
					core::id::AtomID id( j, i );
					//TR << testPose.xyz(id)(1) << " : " << challengePose.xyz(id)(1) << std::endl;
					TS_ASSERT_DELTA( testPose.xyz(id), challengePose.xyz(id), 0.001 );
				}
			}
		}
		TR << "End Test" << std::endl;
	}

	void test_self_graft_recovery_oneFragment_using_end2endMatching()
	{
		using namespace protocols::motif_grafting::movers;

		TR << "Start Test 3" << std::endl;

		//Create the mover
		TR << "Creating mover" << std::endl;
		MotifGraftMover mgm;
		TR << "Intializing mover" << std::endl;

		// initialize the parameters
		std::string s_context     = "protocols/motif_grafting/2nm1_context.pdb";       //context structure
		std::string s_motif      = "protocols/motif_grafting/1zx3a_as_motif_1fragment.pdb"; //motif with one fragment
		core::Real r_RMSD      = 0.5;               //RMSD tol
		core::Real r_NC_RMSD     = 0.5;               //NC points RMSD tol
		core::Size  s_ClashTol    = 5;                //Clash Score tol
		std::string s_comb_d     = "0:0";              //Combinatory delta
		std::string s_max_frag_d   = "0:0";              //Max Fragment relpacement delta
		std::string s_clash_res    = "GLY";              //Clash score test residue
		std::string s_hotspots    = "";                //hotspots
		bool   b_full_bb_alignment  = false;              //full_bb_alignment
		bool   b_allow_ind_frag_alig = true;               //allow_independent_alignment_per_fragment
		bool   b_only_copy_sidechain = false;              //only copy sidechains for hotspots (opposed to graft)
		bool  b_only_if_Npoint_eq  = false;              //graft only if Npoint aminoacid identities are equal
		bool  b_only_if_Cpoint_eq  = false;              //graft only if Cpoint aminoacid identities are equal
		bool   b_revert_to_native  = false;              //revert to the native scafold sequence after grafting
		bool   b_allow_repeat    = false;              //Allow repeat same graft

		mgm.init_parameters(
			s_context,
			s_motif,
			r_RMSD,
			r_NC_RMSD,
			s_ClashTol,
			s_comb_d,
			s_max_frag_d,
			s_clash_res,
			s_hotspots,
			b_full_bb_alignment,
			b_allow_ind_frag_alig,
			b_only_copy_sidechain,
			b_only_if_Npoint_eq,
			b_only_if_Cpoint_eq,
			b_revert_to_native,
			b_allow_repeat);

		//Import the test scaffold
		TR << "Importing Pose" << std::endl;
		core::pose::Pose testPose;
		core::import_pose::pose_from_file( testPose, "protocols/motif_grafting/1zx3a_translated_and_rotated.pdb" , core::import_pose::PDB_file);

		//Get the original sequence
		std::string seq_motif_before=testPose.chain_sequence(1);

		TR << "Applying mover" << std::endl;
		//Apply the mover
		mgm.apply(testPose);

		//get the sequence of the epigraft
		std::string seq_motif_after=testPose.chain_sequence(2);

		//Check that the sequences are equal:
		TS_ASSERT(seq_motif_before.size()==seq_motif_after.size());
		for ( core::Size i=0; i < seq_motif_before.size(); ++i ) {
			//TR << seq_motif_before.at(i) << " : " << seq_motif_after.at(i) << std::endl;
			TS_ASSERT( seq_motif_before.at(i) == seq_motif_after.at(i) );
		}

		/////Import the target scaffold "Real Answer"
		core::pose::Pose challengePose;
		core::import_pose::pose_from_file( challengePose, "protocols/motif_grafting/2nm1_1zx3a_hybrid.pdb" , core::import_pose::PDB_file);


		/////Check that "Answer" and the "Real Answer" sequences are equal size:
		TS_ASSERT( testPose.total_residue() == challengePose.total_residue() );

		/////Check that the atomic positions of the "Answer" and the "Real Answer" are equal
		for ( core::Size i = 1; i <= testPose.total_residue(); ++i ) {
			//core::conformation::ResidueOP resA = new core::conformation::Residue (testPose.residue( i ) );
			//core::conformation::ResidueOP resB = new core::conformation::Residue (challengePose.residue( i ) );
			for ( core::Size j = 1; j <= testPose.residue_type(i).natoms(); ++j ) {
				//This is because of the hack on the 1st hydrogen (1H) of the first motif
				if ( !testPose.residue(i).atom_is_hydrogen(j) ) {
					core::id::AtomID id( j, i );
					//TR << testPose.xyz(id)(1) << " : " << challengePose.xyz(id)(1) << std::endl;
					TS_ASSERT_DELTA( testPose.xyz(id), challengePose.xyz(id), 0.001 );
				}
			}
		}
		TR << "End Test" << std::endl;
	}


	void test_self_graft_recovery_twoFragment_using_end2endMatching()
	{
		using namespace protocols::motif_grafting::movers;

		TR << "Start Test 4" << std::endl;

		//Create the mover
		TR << "Creating mover" << std::endl;
		MotifGraftMover mgm;
		TR << "Intializing mover" << std::endl;

		// initialize the parameters
		std::string s_context     = "protocols/motif_grafting/2nm1_context.pdb";       //context structure
		std::string s_motif      = "protocols/motif_grafting/1zx3a_as_motif_2fragment.pdb"; //motif with one fragment
		core::Real r_RMSD      = 0.5;               //RMSD tol
		core::Real r_NC_RMSD     = 0.5;               //NC points RMSD tol
		core::Size  s_ClashTol    = 5;                //Clash Score tol
		std::string s_comb_d     = "0:0,0:0";            //Combinatory delta
		std::string s_max_frag_d   = "0:0,0:0";            //Max Fragment relpacement delta
		std::string s_clash_res    = "GLY";              //Clash score test residue
		std::string s_hotspots    = "";                //hotspots
		bool   b_full_bb_alignment  = true;               //full_bb_alignment
		bool   b_allow_ind_frag_alig = true;               //allow_independent_alignment_per_fragment
		bool   b_only_copy_sidechain = true;               //only copy sidechains for hotspots (opposed to graft)
		bool  b_only_if_Npoint_eq  = true;               //graft only if Npoint aminoacid identities are equal
		bool  b_only_if_Cpoint_eq  = true;               //graft only if Cpoint aminoacid identities are equal
		bool   b_revert_to_native  = true;               //revert to the native scafold sequence after grafting
		bool   b_allow_repeat    = false;              //Allow repeat same graft

		mgm.init_parameters(
			s_context,
			s_motif,
			r_RMSD,
			r_NC_RMSD,
			s_ClashTol,
			s_comb_d,
			s_max_frag_d,
			s_clash_res,
			s_hotspots,
			b_full_bb_alignment,
			b_allow_ind_frag_alig,
			b_only_copy_sidechain,
			b_only_if_Npoint_eq,
			b_only_if_Cpoint_eq,
			b_revert_to_native,
			b_allow_repeat);

		//Import the test scaffold
		TR << "Importing Pose" << std::endl;
		core::pose::Pose testPose;
		core::import_pose::pose_from_file( testPose, "protocols/motif_grafting/1zx3a_translated_and_rotated.pdb" , core::import_pose::PDB_file);

		//Get the original sequence
		std::string seq_motif_before=testPose.chain_sequence(1);

		TR << "Applying mover" << std::endl;
		//Apply the mover
		mgm.apply(testPose);

		//get the sequence of the epigraft
		std::string seq_motif_after=testPose.chain_sequence(2);

		//Check that the sequences are equal:
		TS_ASSERT(seq_motif_before.size()==seq_motif_after.size());
		for ( core::Size i=0; i < seq_motif_before.size(); ++i ) {
			//TR << seq_motif_before.at(i) << " : " << seq_motif_after.at(i) << std::endl;
			TS_ASSERT( seq_motif_before.at(i) == seq_motif_after.at(i) );
		}

		/////Import the target scaffold "Real Answer"
		core::pose::Pose challengePose;
		core::import_pose::pose_from_file( challengePose, "protocols/motif_grafting/2nm1_1zx3a_hybrid.pdb" , core::import_pose::PDB_file);


		/////Check that "Answer" and the "Real Answer" sequences are equal size:
		TS_ASSERT( testPose.total_residue() == challengePose.total_residue() );

		/////Check that the atomic positions of the "Answer" and the "Real Answer" are equal
		for ( core::Size i = 1; i <= testPose.total_residue(); ++i ) {
			//core::conformation::ResidueOP resA = new core::conformation::Residue (testPose.residue( i ) );
			//core::conformation::ResidueOP resB = new core::conformation::Residue (challengePose.residue( i ) );
			for ( core::Size j = 1; j <= testPose.residue_type(i).natoms(); ++j ) {
				//This is because of the hack on the 1st hydrogen (1H) of the first motif
				if ( !testPose.residue(i).atom_is_hydrogen(j) ) {
					core::id::AtomID id( j, i );
					//TR << testPose.xyz(id)(1) << " : " << challengePose.xyz(id)(1) << std::endl;
					TS_ASSERT_DELTA( testPose.xyz(id), challengePose.xyz(id), 0.001 );
				}
			}
		}
		TR << "End Test" << std::endl;
	}

};
}

