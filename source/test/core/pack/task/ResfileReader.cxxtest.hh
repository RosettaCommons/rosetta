// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/ResfileReader.cxxtest.hh
/// @brief  test suite for resfile reader
/// @author Steven Lewis

// Test headers
#include <cxxtest/TestSuite.h>

#include <core/types.hh>

#include <test/core/init_util.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/RotamerSampleOptions.hh>

#include <core/pack/task/ResfileReader.hh>

#include <string>
#include <sstream> //stringstreams can convert integers into strings type-safely for comparisons en masse

//Auto Headers
#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>

//I'm lazy using's
using namespace core;
using namespace pack;
using namespace task;
using namespace pose;
using namespace chemical;
using std::string;
using std::stringstream;

// --------------- Test Class --------------- //

class ResfileReaderTests : public CxxTest::TestSuite {

public:

	Pose pose;
	PackerTaskOP the_task;
	bool cache_option__interactive_;

	// --------------- Fixtures --------------- //

	void setUp() {
		core_init();
		core::import_pose::pose_from_file( pose, "core/pack/task/resfile_test.pdb" , core::import_pose::PDB_file);
		the_task = TaskFactory::create_packer_task( pose );

		cache_option__interactive_ = basic::options::option[ basic::options::OptionKeys::run::interactive ].value();
		basic::options::option[ basic::options::OptionKeys::run::interactive ].value(true);

	}

	void tearDown() {
		pose.clear(); //nuke that sucker in case it got altered
		//nothing necessary for OP packer task - setUp's regeneration of a fresh copy will cause the old OP's data to delete itself

		basic::options::option[ basic::options::OptionKeys::run::interactive ].value(cache_option__interactive_);

	}

	// ------------- Helper Functions ------------- //

	//extracts the restype list as a string of one-letter name codes
	string decompose_restypes_list(Size resid){
		string oneletters;
		for ( ResidueLevelTask::ResidueTypeCOPListConstIter
				restype_iter = the_task->residue_task( resid ).allowed_residue_types_begin(),
				restype_iter_end = the_task->residue_task( resid ).allowed_residue_types_end();
				restype_iter != restype_iter_end; ++restype_iter ) {

			oneletters += ((**restype_iter).name1());
		}//for

		return oneletters;
	}//decompose_restypes_list

	//extracts the extrachi levels into a string for ease of comparison; format in test function
	string decompose_EX_info( Size resid, ResidueTypeCOP const phe_op, ResidueTypeCOP const ala_op){
		//stringstreams allow coercion of ExtraRotSample -> int -> char -> string for ease of comparison
		//each position in the returned string represents a different query (see below)
		std::ostringstream EXstream;
		ResidueType const & phe = *phe_op;
		ResidueType const & ala = *ala_op;

		EXstream << the_task->residue_task( resid ).extrachi_sample_level(true, 1, phe) //position 0 in string
			<< the_task->residue_task( resid ).extrachi_sample_level(false, 1, phe)//1
			<< the_task->residue_task( resid ).extrachi_sample_level(true, 2, phe) //2
			<< the_task->residue_task( resid ).extrachi_sample_level(false, 2, phe)//3
			<< the_task->residue_task( resid ).extrachi_sample_level(true, 3, phe) //4
			<< the_task->residue_task( resid ).extrachi_sample_level(false, 3, phe)//5
			<< the_task->residue_task( resid ).extrachi_sample_level(true, 4, phe) //6
			<< the_task->residue_task( resid ).extrachi_sample_level(false, 4, phe)//7
			<< the_task->residue_task( resid ).extrachi_sample_level(true, 1, ala) //8
			<< the_task->residue_task( resid ).extrachi_sample_level(false, 1, ala)//9
			<< the_task->residue_task( resid ).extrachi_sample_level(true, 2, ala) //10
			<< the_task->residue_task( resid ).extrachi_sample_level(false, 2, ala)//11
			<< the_task->residue_task( resid ).extrachi_sample_level(true, 3, ala) //12
			<< the_task->residue_task( resid ).extrachi_sample_level(false, 3, ala)//13
			<< the_task->residue_task( resid ).extrachi_sample_level(true, 4, ala) //14
			<< the_task->residue_task( resid ).extrachi_sample_level(false, 4, ala);//15

		return EXstream.str();
	}//decompose_EX_info

	// --------------- Test Cases --------------- //

	//this test will test the "primary modes" ALLAAxc, NOTAA, etc
	void test_resfile_primary_modes() {

		parse_resfile(pose, *the_task, "core/pack/task/resfile_test.resfile" );

		//check that some are being designed
		TS_ASSERT( the_task->design_any());//, true);

		//check designs positively for ALLAAwc, ALLAAxc, PIKAA, NOTAA, POLAR, APOLAR
		TS_ASSERT( the_task->design_residue(1));//, true);
		TS_ASSERT( the_task->design_residue(8));//, true);
		TS_ASSERT( the_task->design_residue(7));//, true);
		TS_ASSERT( the_task->design_residue(4));//, true);
		TS_ASSERT( the_task->design_residue(11));//, true);
		TS_ASSERT( the_task->design_residue(15));//, true);

		//spot-check negatives; one from default set (NATRO and NATAA)
		TS_ASSERT( !the_task->design_residue(5));//), false);
		TS_ASSERT( !the_task->design_residue(17));//), false);
		TS_ASSERT( !the_task->design_residue(12));//), false);

		//check proper number being packed
		TS_ASSERT_EQUALS( the_task->num_to_be_packed(), 14u);

		//check packs positively for ALLAAwc, ALLAAxc, PIKAA, NOTAA, POLAR, APOLAR, NATAA
		TS_ASSERT( the_task->pack_residue(1));//), true);
		TS_ASSERT( the_task->pack_residue(8));//), true);
		TS_ASSERT( the_task->pack_residue(7));//), true);
		TS_ASSERT( the_task->pack_residue(4));//), true);
		TS_ASSERT( the_task->pack_residue(11));//), true);
		TS_ASSERT( the_task->pack_residue(15));//), true);
		TS_ASSERT( the_task->pack_residue(5));//), true);

		//spot-check negatives; one from default set for NATRO
		TS_ASSERT( !the_task->pack_residue(17));//), false);
		TS_ASSERT( !the_task->pack_residue(9));//), false);

		////////////////////////////////////////////////////////////////////////
		//now testing the contents of individual residue records
		//this uses simple string comparisons on purpose
		//altering the sorting in the PackerTask_ WILL cause this to fail
		//this is deliberate, because it will require to you recheck the sorting in
		//the PackerTask_ if you change it, and then input that here
		//SML

		//uncomment the for loop to see the contents of each residue types list
		/////////////////////////////////////////////////////////////////////////


		//   for ( int q = 1; q < 30 /*assuming pose size is 30 */; q++ ) {
		//    std::cout << q << "  " << decompose_restypes_list(q) << std::endl;
		//   }


		TS_ASSERT_EQUALS(decompose_restypes_list(1), "ACDEFGHHIKLMNPQRSTVWY");  //ALLAAwc
		TS_ASSERT_EQUALS(decompose_restypes_list(2), "E"); //PIKAA E
		TS_ASSERT_EQUALS(decompose_restypes_list(3), ""); //NATRO means empty string (no rotamers)
		TS_ASSERT_EQUALS(decompose_restypes_list(4), "ACDEFGHHKLMNP"); //NOTAA QRSTVWYI
		TS_ASSERT_EQUALS(decompose_restypes_list(5), "I"); //NATAA when natively I
		TS_ASSERT_EQUALS(decompose_restypes_list(7), "ACDEF"); //PIKAA ACDEF
		TS_ASSERT_EQUALS(decompose_restypes_list(8), "ADEFGHHIKLMNPQRSTVWY"); //ALLAAxc
		TS_ASSERT_EQUALS(decompose_restypes_list(11), "DEHHKNQRST"); //POLAR
		TS_ASSERT_EQUALS(decompose_restypes_list(15), "ACFGILMPVWY"); //APOLAR
		TS_ASSERT_EQUALS(decompose_restypes_list(17), ""); //manual NATRO
		TS_ASSERT_EQUALS(decompose_restypes_list(21), "ADHH"); //PIKAA ADH, for a PIKAA with H
		TS_ASSERT_EQUALS(decompose_restypes_list(22), "ACDEGHHILMNPQRSWY"); //PIKAA QWERYIPLHGDSACNM #out of order
		TS_ASSERT_EQUALS(decompose_restypes_list(23), "ADFGIKMPSY"); //NOTAA NVCLHTREWQ #also out of order
		TS_ASSERT_EQUALS(decompose_restypes_list(24), "ACDGILMNPQSVW"); //PIKAA QWIPGLDSACVNM #histidine removed
		TS_ASSERT_EQUALS(decompose_restypes_list(25), "ADFGHHIKMPSY"); //NOTAA NVCLTREWQ #H removed

		TS_ASSERT(the_task->include_current(26)); //26 use_input_sc
		TS_ASSERT(!the_task->include_current(27));//27 #default

	}// end test_resfile_primary_modes

	//this test will examine the EX commands
	void test_resfile_EX_commands(){

		core::Size extrachi_cutoff_limit(core::pack::task::EXTRACHI_CUTOFF_LIMIT); //we want to cast the int to an unsigned

		parse_resfile(pose, *the_task, "core/pack/task/EX_test.resfile" );

		//we need instances of an aromatic and nonaromatic residue
		//fastest route is filching ResidueTypeSet reference out of pose
		ResidueTypeCOP phe = pose.residue_type(1).residue_type_set()->name_map("PHE").get_self_ptr();
		ResidueTypeCOP ala = pose.residue_type(1).residue_type_set()->name_map("ALA").get_self_ptr();

		//the "magic number" strings represent the expected data
		//the format, presented as (position indexed from 0, meaning) is as follows:
		//0  buried aromatic EX 1 level
		//1  nonburied aromatic EX 1 level
		//2  buried aromatic EX 2 level
		//3  nonburied aromatic EX 2 level
		//4  buried EX 3 level
		//5  nonburied EX 3 level
		//6  buried EX 4 level
		//7  nonburied EX 4 level
		//8  buried nonaromatic EX 1 level
		//9  nonburied nonaromatic EX 1 level
		//10 buried nonaromatic EX 2 level
		//11 nonburied nonaromatic EX 2 level
		//12  buried nonaromatic EX 3 level
		//13  nonburied nonaromatic EX 3 level
		//14 buried nonaromatic EX 4 level
		//15 nonburied nonaromatic EX 4 level
		TS_ASSERT_EQUALS(decompose_EX_info(1, phe, ala),  "1000000010000000"); //1 EX 1
		TS_ASSERT_EQUALS(decompose_EX_info(2, phe, ala),  "1010000010100000"); //2 EX 1 EX 2
		TS_ASSERT_EQUALS(decompose_EX_info(3, phe, ala),  "1010000010100000"); //3 EX 2 EX 1
		TS_ASSERT_EQUALS(decompose_EX_info(4, phe, ala),  "1010100010101000"); //4 EX 1 EX 2 EX 3
		TS_ASSERT_EQUALS(decompose_EX_info(5, phe, ala),  "0000101000001010"); //5 EX 3 EX 4
		TS_ASSERT_EQUALS(decompose_EX_info(6, phe, ala),  "1010101010101010"); //6 EX 1 EX 2 EX 4 EX 3
		TS_ASSERT_EQUALS(decompose_EX_info(7, phe, ala),  "1000000000000000"); //7 EX ARO 1
		TS_ASSERT_EQUALS(decompose_EX_info(8, phe, ala),  "0010000000000000"); //8 EX ARO 2
		TS_ASSERT_EQUALS(decompose_EX_info(9, phe, ala),  "1010000000000000"); //9 EX ARO 1 EX ARO 2
		TS_ASSERT_EQUALS(decompose_EX_info(10, phe, ala), "1010001000100010"); //10 EX ARO 1 EX 2 EX 4
		TS_ASSERT_EQUALS(decompose_EX_info(11, phe, ala), "4000000040000000"); //11 EX 1 LEVEL 4
		TS_ASSERT_EQUALS(decompose_EX_info(12, phe, ala), "0020000000200000"); //12 EX 2 LEVEL 2
		TS_ASSERT_EQUALS(decompose_EX_info(13, phe, ala), "0000500000005000"); //13 EX 3 LEVEL 5
		TS_ASSERT_EQUALS(decompose_EX_info(14, phe, ala), "3040201000102010"); //14 EX ARO 1 LEVEL 3 EX 3 LEVEL 2 EX 2 EX ARO 2 LEVEL 4 EX 4
		TS_ASSERT_EQUALS(decompose_EX_info(15, phe, ala), "3000000030000000"); //15 EX 1 LEVEL 2 EX 1 LEVEL 3
		TS_ASSERT_EQUALS(decompose_EX_info(16, phe, ala), "2000000010000000"); //16 EX ARO 1 LEVEL 2 EX 1
		TS_ASSERT_EQUALS(decompose_EX_info(17, phe, ala), "0000006000000060"); //17 EX 4 LEVEL 6
		TS_ASSERT_EQUALS(decompose_EX_info(18, phe, ala), "1010000010100000"); //18 is default
		TS_ASSERT_EQUALS(decompose_EX_info(31, phe, ala), "1000000010000000"); //this is testing that resid 31 tracks PDB id 51 _

		TS_ASSERT_EQUALS(the_task->residue_task(24).extrachi_cutoff(), extrachi_cutoff_limit); //24 is default
		TS_ASSERT_EQUALS(the_task->residue_task(25).extrachi_cutoff(), 3u); //25 EX_cutoff 3
		TS_ASSERT_EQUALS(the_task->residue_task(26).extrachi_cutoff(), 0u); //26 EX_cutoff 0
		TS_ASSERT_EQUALS(the_task->residue_task(27).extrachi_cutoff(), extrachi_cutoff_limit); //27 EX_cutoff 99999999 #should leave at 18 = EXTRACHI_CUTOFF_LEVEL


	}//end test_resfile_EX_commands

	// void  test_generate_output_resfiles(){
	//  std::cout << "generate_output_resfiles" << std::endl;
	//  //std::string directory( "/Users/momeara/mini_pure2/test/core/pack/task/test_resfiles/"); //"core/pack/task/test_resfiles");
	//  std::string directory( "core/pack/task/test_resfiles/");
	//  utility::vector1<std::string> filenames;
	//  if ( utility::file::list_dir( directory, filenames ) ){
	//   std::cout << "error reading directory " << directory << "." << std::endl;
	//   TS_ASSERT(false);
	//  }
	//  for( utility::vector1<std::string>::iterator input_filename = filenames.begin(); input_filename != filenames.end(); ++input_filename ){
	//   if ( utility::file::file_extension( *input_filename ) != "resfile") continue;
	//   std::string output_filename( utility::file::file_basename( *input_filename ) + ".output" );
	//   std::cout << "generating output '" << directory << output_filename << "' from input file '" << directory << *input_filename << "'." << std::endl;
	//   the_task = TaskFactory::create_packer_task( pose );
	//   the_task->read_resfile( directory + *input_filename );
	//
	//   std::string actual_output = the_task->task_string( pose );
	//   utility::file::file_delete( directory + output_filename);
	//   utility::file::create_blank_file( directory + output_filename );
	//   std::ofstream output_file;
	//   output_file.open( (directory+output_filename).c_str() );
	//   output_file.close();
	//  }
	// }
	//
	// void generate_output_resfile( std::string const & directory,
	//                std::string const & input_filename ){
	//
	//  if ( utility::file::file_extension( input_filename ) != "resfile"){
	//   std::cout << "The output should have the extension 'resfile'" << std::endl;
	//   TS_ASSERT(false);
	//  }
	//
	//  std::string output_filename( utility::file::file_basename( input_filename ) + ".output" );
	//  std::cout << "generating output '" << directory << output_filename << "' from input file '" << directory << input_filename << "'." << std::endl;
	//  the_task = TaskFactory::create_packer_task( pose );
	//  the_task->read_resfile( directory + input_filename );
	//
	//  std::string actual_output = the_task->task_string( pose );
	//  utility::file::file_delete( directory + output_filename);
	//  utility::file::create_blank_file( directory + output_filename );
	//  std::ofstream output_file;
	//  output_file.open( (directory+output_filename).c_str() );
	//  output_file  << actual_output;
	//  output_file.close();
	// }
	//
	// void compare_input_output( std::string const & directory,
	//            std::string const & input_filename ){
	//
	//  std::string output_filename( utility::file::file_basename( input_filename ) + ".output" );
	//  the_task = TaskFactory::create_packer_task( pose );
	//  the_task->read_resfile( directory + input_filename );
	//  std::string actual_output = the_task->task_string( pose );
	//  std::string expected_output;
	//  utility::io::izstream file( directory + output_filename );
	//  if (!file) {
	//   std::cout << "Cannot open file " << output_filename << " run generate_output_resfile in ResfileReader.cxxtest.hh to generate";
	//   TS_ASSERT(false);
	//  }
	//  utility::slurp( file, expected_output );
	//  TS_ASSERT_EQUALS( actual_output, expected_output );
	// }
	//
	// void test_output_task_string(){
	//  std::string directory( "core/pack/task/test_resfiles/"); //"core/pack/task/test_resfiles");
	//  std::cout << "test_output_task_string, putting output in directory " << directory << std::endl;
	//
	//  //for resfile in directory, read in resfile, compare it with expected output
	//  utility::vector1<std::string> filenames;
	//  if ( utility::file::list_dir( directory, filenames ) ) {
	//   std::cout << "error reading directory " << directory << "." << std::endl;
	//   TS_ASSERT(false);
	//  }
	//
	//  for( utility::vector1<std::string>::iterator input_it=filenames.begin(); input_it !=filenames.end(); ++input_it ){
	//   if ( utility::file::file_extension( *input_it ) != "resfile") continue;
	//   generate_output_resfile( directory, *input_it );
	//   compare_input_output( directory, *input_it );
	//  }
	// }


	void test_resid_parsing(){
		{
			PackerTaskOP ptask = TaskFactory::create_packer_task(pose);
			stringstream resfile; resfile << "NATRO\nSTART\n1 - 3 _ ALLAA";
			parse_resfile_string(pose, *ptask, "dummy_filename", resfile.str() );
			TS_ASSERT(ptask->pack_residue(1));
			TS_ASSERT(ptask->pack_residue(2));
			TS_ASSERT(ptask->pack_residue(3));
			TS_ASSERT(!ptask->pack_residue(4));
			TS_ASSERT(!ptask->pack_residue(5));
		}
		{
			PackerTaskOP ptask = TaskFactory::create_packer_task(pose);
			stringstream resfile; resfile << "NATRO\nSTART\n2- 4  _ ALLAA";
			try{
				parse_resfile_string(pose, *ptask, "dummy_filename", resfile.str() );
				TS_FAIL("Didn't catch malformed range resid.");
			} catch(ResfileReaderException) {}
		}
		{
			PackerTaskOP ptask = TaskFactory::create_packer_task(pose);
			stringstream resfile; resfile << "NATRO\nSTART\n3 - 5 _ ALLAA";
			parse_resfile_string(pose, *ptask, "dummy_filename", resfile.str() );
			TS_ASSERT(!ptask->pack_residue(1));
			TS_ASSERT(!ptask->pack_residue(2));
			TS_ASSERT(ptask->pack_residue(3));
			TS_ASSERT(ptask->pack_residue(4));
			TS_ASSERT(ptask->pack_residue(5));
		}
		{
			PackerTaskOP ptask = TaskFactory::create_packer_task(pose);
			stringstream resfile; resfile << "NATRO\nSTART\n* _ ALLAA";
			parse_resfile_string(pose, *ptask, "dummy_filename", resfile.str() );
			for ( Size i = 1; i <= pose.total_residue(); ++i ) {
				TS_ASSERT(ptask->pack_residue(i));
			}
		}
		{
			PackerTaskOP ptask = TaskFactory::create_packer_task(pose);
			stringstream resfile; resfile << "START\n* _ NATRO\n1 - 3 _ NATAA";
			parse_resfile_string(pose, *ptask, "dummy_filename", resfile.str() );
			TS_ASSERT(ptask->pack_residue(1));
			TS_ASSERT(ptask->pack_residue(2));
			TS_ASSERT(ptask->pack_residue(3));
			for ( Size i = 4; i <= pose.total_residue(); ++i ) {
				TS_ASSERT(!ptask->pack_residue(i));
			}
		}
		{
			PackerTaskOP ptask = TaskFactory::create_packer_task(pose);
			stringstream resfile; resfile << "START\n1 - 3 _ NATAA\n* _ NATRO";
			parse_resfile_string(pose, *ptask, "dummy_filename", resfile.str() );
			TS_ASSERT(ptask->pack_residue(1));
			TS_ASSERT(ptask->pack_residue(2));
			TS_ASSERT(ptask->pack_residue(3));
			for ( Size i = 4; i <= pose.total_residue(); ++i ) {
				TS_ASSERT(!ptask->pack_residue(i));
			}
		}
		{
			stringstream resfile; resfile << "NATRO\nSTART\n * - 3 _ NATRO";
			try{
				ResfileContents(pose, "dummy_filename", resfile);
				TS_FAIL("Didn't catch bad chain and range resid.");
			} catch(ResfileReaderException){
			}
		}
		{
			stringstream resfile; resfile << "START\n2 _ ";
			try{
				ResfileContents(pose, "dummy_filename", resfile);
				TS_FAIL("Didn't catch bad missing commands for single resid.");
			} catch(ResfileReaderException){}
		}
		{
			stringstream resfile; resfile << "START\n 3 - 4 _ ";
			try{
				ResfileContents(pose, "dummy_filename", resfile);
				TS_FAIL("Didn't catch bad missing commands for range resid.");
			} catch(ResfileReaderException){}
		}
		{
			stringstream resfile; resfile << "START\n * _ ";
			try{
				ResfileContents(pose, "dummy_filename", resfile);
				TS_FAIL("Didn't catch bad missing commands for chain resid.");
			} catch(ResfileReaderException){}
		}
		{
			stringstream resfile; resfile << "START\n3 - 2  _ NATAA";
			try{
				ResfileContents(pose, "dummy_filename", resfile);
				TS_FAIL("Didn't catch that the start residue must come before the end residue in a range resid.");
			} catch(ResfileReaderException){}
		}
		{
			stringstream resfile; resfile << "START\n2 % NATAA";
			try{
				ResfileContents(pose, "dummy_fname", resfile);
				TS_FAIL("Didn't catch that the chain identifier is not in [_A-Za-z].");
			} catch(ResfileReaderException){}
		}
		{
			stringstream resfile; resfile << "START\n2 __ NATAA";
			try{
				ResfileContents(pose, "dummy_filename", resfile);
				TS_FAIL("Didn't catch that the chain must be just a single character.");
			} catch(ResfileReaderException){}
		}
	}

};//end class
