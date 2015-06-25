// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/reference_pose/ReferencePose.cxxtest.hh
/// @brief  Unit tests for the ReferencePose and ReferencePoseSet classes.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>

#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/reference_pose/ReferencePoseSet.fwd.hh>
#include <core/pose/reference_pose/ReferencePose.fwd.hh>
#include <core/pose/reference_pose/ReferencePoseSet.hh>
#include <core/pose/reference_pose/ReferencePose.hh>
#include <utility/vector1.hh>

#include <test/core/init_util.hh> //necessary if there is tracer output
#include <basic/Tracer.hh>

//Auto Headers

static basic::Tracer TR("core.pose.reference_pose.ReferencePose.cxxtest");

// --------------- Test Class --------------- //

class ReferencePoseTests : public CxxTest::TestSuite {

public:

	// typedefs
	typedef core::pose::reference_pose::ReferencePose ReferencePose;
	typedef core::pose::reference_pose::ReferencePoseSet ReferencePoseSet;

	// shared initialization
	void setUp() {
		core_init(); //necessary if there is tracer output for a failure
	}

	// shared finalization
	void tearDown() {
	}

	// ------------- Helper Functions ------------- //


	// --------------- Test Cases --------------- //

	void test_ReferencePoseSet_setup() {
		TR << "Starting ReferencePoseSet_setup()." << std::endl;
		TR << "This test is intended to test whether a ReferencePose can be created properly." << std::endl;
		TR << "Failure means that we were unable to create the ReferencePose object in the" << std::endl;
		TR << "pose object, or that refence pose numbers were not properly initialized." << std::endl;
		TR << "Importing 1AFO_AB.pdb." << std::endl;
		core::pose::PoseOP pose(new core::pose::Pose);
		core::import_pose::pose_from_pdb( *pose, "core/pose/reference_pose/1AFO_AB.pdb" );
		
		TR << "Setting up reference pose from imported pose." << std::endl;
		pose->reference_pose_from_current("refpose1");
		
		TR << "Checking ReferencePose inital state. We expect there to be an entry for each residue, mapping that residue index onto itself." << std::endl;
		TR << "Residue\tMap value" << std::endl;
		for(core::Size ir=1, irmax=pose->n_residue(); ir<=irmax; ++ir) { //Loop through all residues in the pose.
			core::Size const cor_ir( pose->corresponding_residue_in_current( ir, "refpose1" ) );
			TR << ir << "\t" << cor_ir << std::endl;
			TS_ASSERT( cor_ir==ir );
		}
		
		TR.flush();
		return;
	}
	
	void test_ReferencePoseSet_deleting() {
		TR << "Starting ReferencePoseSet_deleting()." << std::endl;
		TR << "This test is intended to test whether a ReferencePose properly maintains residue" << std::endl;
		TR << "indices as residues are deleted.  Failure means that we have not properly shifted" << std::endl;
		TR << "residue indices when we deleted a residue." << std::endl;
		TR << "Importing 1AFO_AB.pdb." << std::endl;
		core::pose::PoseOP pose(new core::pose::Pose);
		core::import_pose::pose_from_pdb( *pose, "core/pose/reference_pose/1AFO_AB.pdb" );
		
		core::Size const nres_original( pose->n_residue() );
		TR << "There are " << nres_original << " residues in the original pose." << std::endl;
		
		TR << "Setting up reference pose from imported pose." << std::endl;
		pose->reference_pose_from_current("refpose2");
		
		TR << "Deleting residue 7." << std::endl;
		pose->delete_polymer_residue(7);
		
		TR << "Checking ReferencePose state after deletion. We expect there to be an entry for each residue in the original pose, mapping that residue index onto itself for residues prior to 7, onto 0 for residue 7, and onto residue n-1 for residues past 7." << std::endl;
		TR << "Residue\tMap value" << std::endl;
		for(core::Size ir=1; ir<=nres_original; ++ir) { //Loop through all residues in the pose.
			core::Size const cor_ir( pose->corresponding_residue_in_current( ir, "refpose2" ) );
			TR << ir << "\t" << cor_ir << std::endl;
			if(ir<7) {
				TS_ASSERT( cor_ir==ir );
			}
			else if (ir==7) {
				TS_ASSERT( cor_ir==0 );
			}
			else if (ir>7) {
				TS_ASSERT( cor_ir==ir-1 );
			}
		}
		
		TR.flush();
		return;
	}
	
	void test_ReferencePose_string_parsing() {
		TR << "Starting ReferencePose_string_parsing()." << std::endl;
		TR << "This test tests whether a string referring to a stored ReferencePose can be parsed" << std::endl;
		TR << "properly.  Such strings must be of the form \"refpose(<refpose_name>,<residue_number>)\"" << std::endl;
		TR << "or \"refpose(<refpose_name>,<residue_number>)+/-<offset_value>\"." << std::endl;

		utility::vector1<std::string> str;
		utility::vector1<bool> result;
		utility::vector1<std::string> result_refpose_name;
		result_refpose_name.resize(6,"");
		utility::vector1<core::Size> result_resnum;
		utility::vector1<signed long> result_offset;
		result_resnum.resize(6,0);
		result_offset.resize(6,0);

		str.push_back( "Bad_input." ); //Bad string
		str.push_back( "refpose(badinput" ); //Bad string
		str.push_back( "refpose(refpose1,)" ); //Bad string
		str.push_back( "refpose(refpose2,9)" ); //Good string
		str.push_back( "refpose(,10)" ); //Bad string
		str.push_back( "refpose(refpose3,11)-3" ); //Good string
		
		TR << "Input_string\tGood_string?\tRefpose_name\tRef_number\tRef_offset" << std::endl;
		for(core::Size i=1; i<=6; ++i)
		{
			result.push_back( core::pose::is_referencepose_number( str[i], result_refpose_name[i], result_resnum[i], result_offset[i] ) );
			TR << str[i] << "\t" << (result[i]?"true":"false") << "\t" << (result_refpose_name[i]==""?"NONE":result_refpose_name[i]) << "\t" << result_resnum[i] << "\t" << result_offset[i] << std::endl;
		}
		
		//Check that the expected failures occurred:
		TS_ASSERT( result[1]==false );
		TS_ASSERT( result[2]==false );
		TS_ASSERT( result[3]==false );
		TS_ASSERT( result[4]==true );
		TS_ASSERT( result[5]==false );
		TS_ASSERT( result[6]==true );

		//Check that the names were extracted properly:
		TS_ASSERT( result_refpose_name[4]=="refpose2" );
		TS_ASSERT( result_refpose_name[6]=="refpose3" );

		//Check that the numbers were extracted properly:
		TS_ASSERT( result_resnum[4]==9 );
		TS_ASSERT( result_resnum[6]==11 );
		TS_ASSERT( result_offset[4]==0 );
		TS_ASSERT( result_offset[6]==-3 );
		
		TR.flush();
		return;
	}

};//end class
