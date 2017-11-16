// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    test/core/io/pdb/file_data.cxxtest.hh
/// @brief   Test suite for StructFileRep methods and utility functions
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <core/io/pdb/Field.hh>
#include <core/io/pdb/pdb_reader.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/StructFileRep.hh>
#include <core/io/pose_to_sfr/PoseToStructFileRepConverter.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/import_pose/import_pose.hh>

// Utility headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR("core.io.StructFileRep.cxxtest");

class StructFileRepTests : public CxxTest::TestSuite {
public:  // Standard methods //////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
	}

	// Destruction
	void tearDown()
	{}


public:  // Tests /////////////////////////////////////////////////////////////
	// Input Methods //////////////////////////////////////////////////////////


	/// @brief Tests PDB import when the -extra_res_fa flag is used to specify a noncanonical residue.
	void test_extra_res_fa_flag() {
		// This has to be the first one called because we don't re-initialize core and we really need to get this extra res
		// in the RTS...
		core_init_with_additional_options( "-obey_ENDMDL -constraints_from_link_records -cst_weight 1 -extra_res_fa core/io/pdb/test.params");
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/io/pdb/extra_res_pose.pdb" , core::import_pose::PDB_file);
		TS_ASSERT_EQUALS( pose.size(), 181 );
	}


	// Output Methods /////////////////////////////////////////////////////////
	// Confirm that LinkInformation and SSBondInformation data are created properly from a Pose.
	void test_get_connectivity_annotation_info()
	{
		using namespace core::io;
		using namespace core::io::pdb;
		core_init_with_additional_options( "-obey_ENDMDL -constraints_from_link_records -cst_weight 1");

		TR <<  "Testing get_connectivity_annotation_info() method of StructFileRep."  << std::endl;

		core::pose::Pose pose;

		// 1BH4 is circulin A, a cyclic peptide with three disulfides.
		//core::import_pose::pose_from_file( pose, "core/io/pdb/1BH4.pdb" , core::import_pose::PDB_file);

		// 4TTL is a cyclic peptide with two disulfides.
		core::import_pose::pose_from_file( pose, "core/io/pdb/4TTL.pdb" , core::import_pose::PDB_file);

		//core::import_pose::pose_from_file( pose, "core/io/pdb/1e68_link.pdb" , core::import_pose::PDB_file);

		core::io::pose_to_sfr::PoseToStructFileRepConverter pose_to_sfr;

		TS_ASSERT_EQUALS( pose_to_sfr.sfr()->link_map().size(),   0 );
		pose_to_sfr.get_connectivity_annotation_info( pose );

		TS_ASSERT_EQUALS( pose_to_sfr.sfr()->ssbond_map().size(), 2 );
		TS_ASSERT_EQUALS( pose_to_sfr.sfr()->link_map().size(),   1 );

		// The following lines came directly from 1BH4:
		//SSBOND   1 CYS A    1    CYS A   17                          1555   1555  2.02
		//SSBOND   2 CYS A    5    CYS A   19                          1555   1555  2.02
		//SSBOND   3 CYS A   10    CYS A   24                          1555   1555  2.02
		//LINK         N   CYS A   1                 C   PRO A  30     1555   1555  1.31

		// The following lines came directly from 4TTL:
		//SSBOND   1 CYS A    2    CYS A    8                          1555   1555  2.05
		//SSBOND   2 CYS A    3    CYS A   16                          1555   1555  2.03
		//LINK         N   GLY A   1                 C   GLY A  22     1555   1555  1.34

		// The following lines come directly from 1e68_link.pdb:
		//LINK         N   MET A   1                 C   TRP A  70     1555   1555  1.33

		core::pose::Pose pose2;
		core::import_pose::pose_from_file( pose2, "core/io/pdb/1e68_link.pdb" , core::import_pose::PDB_file);
		core::io::pose_to_sfr::PoseToStructFileRepConverter pose_to_sfr2;

		TS_ASSERT_EQUALS( pose_to_sfr2.sfr()->link_map().size(),   0 );
		pose_to_sfr2.get_connectivity_annotation_info( pose2 );

		TS_ASSERT_EQUALS( pose_to_sfr2.sfr()->ssbond_map().size(), 0 );
		TS_ASSERT_EQUALS( pose_to_sfr2.sfr()->link_map().size(),   1 );
	}

	///@author Steven Lewis smlewi@gmail.com
	//maybe should be in PoseToStructFileRepConverter.cxxtest.hh
	void test_generate_secondary_structure_informations(){

		core::pose::Pose pose;
		//this string will be both our pose and its secstruct string (ensuring they're the same length...)
		std::string const pose_AND_ss = "LLHHHHHHLLHHHHHHLLEEEEEELLLLEEEEEEELEEEEEELEEELELHHHLHHLHHHLHHHHHHHHLLEEEHHLL";
		core::pose::make_pose_from_sequence(pose, pose_AND_ss, "fa_standard");

		//this is dumb but pose has no string function for this
		for ( core::Size i(1); i<pose.size(); ++i ) {
			pose.set_secstruct(i, pose_AND_ss[i-1]);
		}

		///////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////
		//check that NOTHING happens when option is false
		core::io::pose_to_sfr::PoseToStructFileRepConverter pose_to_sfr_no;
		core::io::StructFileRepOptions options_no;

		options_no.set_output_secondary_structure(false);

		//verify no SHEET or HELIX to start with
		TS_ASSERT_EQUALS( pose_to_sfr_no.sfr()->HELIXInformations().size(),   0 );
		TS_ASSERT_EQUALS( pose_to_sfr_no.sfr()->SHEETInformations().size(),   0 );

		//init SFR
		pose_to_sfr_no.init_from_pose(pose, options_no);

		//verify no SHEET or HELIX still, because init didn't add them b/c option
		TS_ASSERT_EQUALS( pose_to_sfr_no.sfr()->HELIXInformations().size(),   0 );
		TS_ASSERT_EQUALS( pose_to_sfr_no.sfr()->SHEETInformations().size(),   0 );

		//////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////
		//Check that SOMETHING happens when option is true
		//////////////////////////////////////////////////////////////////////////////

		core::io::pose_to_sfr::PoseToStructFileRepConverter pose_to_sfr;
		core::io::StructFileRepOptions options_yes;

		options_yes.set_output_secondary_structure(true);
		options_yes.set_do_not_autoassign_SS(true); //need to override this to dodge dssp on our garbage pose

		//verify no SHEET or HELIX to start with
		TS_ASSERT_EQUALS( pose_to_sfr.sfr()->HELIXInformations().size(),   0 );
		TS_ASSERT_EQUALS( pose_to_sfr.sfr()->SHEETInformations().size(),   0 );

		//init SFR
		pose_to_sfr.init_from_pose(pose, options_yes);

		TS_ASSERT_DIFFERS( pose_to_sfr.sfr()->HELIXInformations().size(),   0 );
		TS_ASSERT_DIFFERS( pose_to_sfr.sfr()->SHEETInformations().size(),   0 );

		//LLHHHHHHLLHHHHHHLLEEEEEELLLLEEEEEEELEEEEEELEEELELHHHLHHLHHHLHHHHHHHHLLEEEHHLL
		TS_ASSERT_EQUALS( pose_to_sfr.sfr()->HELIXInformations().size(),   7 );
		TS_ASSERT_EQUALS( pose_to_sfr.sfr()->SHEETInformations().size(),   6 );

		//////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////
		//pick one sheet and one helix to test specifically; rest we will test with string comparisons
		//type is utility::vector1< HELIXInformation >
		auto const helixvec(pose_to_sfr.sfr()->HELIXInformations());
		//helix 2 is 11-16
		core::io::HELIXInformation const & test_helix(helixvec[2]);
		TS_ASSERT_EQUALS(test_helix.helixID, 2);
		TS_ASSERT_EQUALS(test_helix.helix_name, "  2");
		TS_ASSERT_EQUALS(test_helix.name3_1, "HIS");
		TS_ASSERT_EQUALS(test_helix.chainID1, 'A');
		TS_ASSERT_EQUALS(test_helix.seqNum1, 11);
		TS_ASSERT_EQUALS(test_helix.icode1, ' ');
		TS_ASSERT_EQUALS(test_helix.name3_2, "HIS");
		TS_ASSERT_EQUALS(test_helix.chainID2, 'A');
		TS_ASSERT_EQUALS(test_helix.seqNum2, 16);
		TS_ASSERT_EQUALS(test_helix.icode2, ' ');
		TS_ASSERT_EQUALS(test_helix.helixClass, 1);
		TS_ASSERT_EQUALS(test_helix.comment, "");
		TS_ASSERT_EQUALS(test_helix.length, 6);

		//type is utility::vector1< SHEETInformation >
		auto const sheetvec(pose_to_sfr.sfr()->SHEETInformations());
		//LLHHHHHHLLHHHHHHLLEEEEEELLLLEEEEEEELEEEEEELEEELELHHHLHHLHHHLHHHHHHHHLLEEEHHLL
		//sheet 4 is 44-46
		core::io::SHEETInformation const & test_sheet(sheetvec[4]);
		TS_ASSERT_EQUALS(test_sheet.strand_num, 1);
		TS_ASSERT_EQUALS(test_sheet.sheetID, "  4");
		TS_ASSERT_EQUALS(test_sheet.num_strands, 1);
		TS_ASSERT_EQUALS(test_sheet.name3_1, "GLU");
		TS_ASSERT_EQUALS(test_sheet.chainID1, 'A');
		TS_ASSERT_EQUALS(test_sheet.seqNum1, 44);
		TS_ASSERT_EQUALS(test_sheet.icode1, ' ');
		TS_ASSERT_EQUALS(test_sheet.name3_2, "GLU");
		TS_ASSERT_EQUALS(test_sheet.chainID2, 'A');
		TS_ASSERT_EQUALS(test_sheet.seqNum2, 46);
		TS_ASSERT_EQUALS(test_sheet.icode2, ' ');
		TS_ASSERT_EQUALS(test_sheet.strandClass, 0);

		//seems unlikely that this is the best place to test this, but.... /shrug
		//////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////
		//check PDB style output is correct, spacing, etc.
		//see also http://www.wwpdb.org/documentation/file-format-content/format23/sect5.html

		//this is kind of dumb, need StructFileRepOptionsCOP for dump_pdb
		//could make deep copy of options_yes, but we'll just test the option system while we're at it!
		core_init_with_additional_options( "-out:file:output_secondary_structure -out:file:do_not_autoassign_SS");
		core::io::StructFileRepOptionsOP options_op(new core::io::StructFileRepOptions);
		//it was initialize from options system in its ctor

		/*notice it isn't padding the records in this copy b/c the editor strips them,
		these are actually padded to 80 chars
		HELIX    1   1 HIS A    3  HIS A    8  1                                   6
		HELIX    2   2 HIS A   11  HIS A   16  1                                   6
		HELIX    3   3 HIS A   50  HIS A   52  1                                   3
		HELIX    4   4 HIS A   54  HIS A   55  1                                   2
		HELIX    5   5 HIS A   57  HIS A   59  1                                   3
		HELIX    6   6 HIS A   61  HIS A   68  1                                   8
		HELIX    7   7 HIS A   74  HIS A   75  1                                   2
		SHEET    1   1 1 GLU A  19  GLU A  24  0
		SHEET    1   2 1 GLU A  29  GLU A  35  0
		SHEET    1   3 1 GLU A  37  GLU A  42  0
		SHEET    1   4 1 GLU A  44  GLU A  46  0
		SHEET    1   5 1 GLU A  48  GLU A  48  0
		SHEET    1   6 1 GLU A  71  GLU A  73  0
		*/

		std::string const right_answer =
			"HELIX    1   1 HIS A    3  HIS A    8  1                                   6    \nHELIX    2   2 HIS A   11  HIS A   16  1                                   6    \nHELIX    3   3 HIS A   50  HIS A   52  1                                   3    \nHELIX    4   4 HIS A   54  HIS A   55  1                                   2    \nHELIX    5   5 HIS A   57  HIS A   59  1                                   3    \nHELIX    6   6 HIS A   61  HIS A   68  1                                   8    \nHELIX    7   7 HIS A   74  HIS A   75  1                                   2    \nSHEET    1   1 1 GLU A  19  GLU A  24  0                                        \nSHEET    1   2 1 GLU A  29  GLU A  35  0                                        \nSHEET    1   3 1 GLU A  37  GLU A  42  0                                        \nSHEET    1   4 1 GLU A  44  GLU A  46  0                                        \nSHEET    1   5 1 GLU A  48  GLU A  48  0                                        \nSHEET    1   6 1 GLU A  71  GLU A  73  0                                        \n";
		std::ostringstream pdb_stream;
		core::io::pdb::dump_pdb( pose, pdb_stream, options_op);
		//std::cout << pdb_stream.str() << std::endl;

		std::istringstream pdb_result(pdb_stream.str());
		std::string test_answer("");
		//convert the first 13 lines of the PDB into our string to compare to
		for ( core::Size lineNum(1); lineNum <= 13 /*13 SS elements*/; ++lineNum ) {
			std::string templine;
			std::getline(pdb_result, templine);
			test_answer += templine + "\n";
		}

		TS_ASSERT_EQUALS(test_answer, right_answer);
	}


	///@author Steven Lewis smlewi@gmail.com
	//maybe should be in PoseToStructFileRepConverter.cxxtest.hh
	//original implementation of this system sorted HELIX and SHEET records in a way that didn't match the ATOM records; this wasn't exposed by the test above.  This function constructs a larger multichain pose to test its improved sorting behavior.
	void test_generate_secondary_structure_informations_twochains(){

		//pose 1
		core::pose::Pose pose1;
		//this string will be both our pose and its secstruct string (ensuring they're the same length...)
		std::string const pose1_AND_ss = "LLHHHHHHLLHHHHHHLLEEEEEELLLLEEEEEEELEEEEEELEEELELHHHLHHLHHHLHHHHHHHHLLEEEHHLL";
		core::pose::make_pose_from_sequence(pose1, pose1_AND_ss, "fa_standard");

		core::pose::Pose pose2;
		std::string const pose2_AND_ss("LLHEEELLLLLLEEEEEELLLLLHHHHLLLLHHHHLLLEEELLEEEELLLEEEELLLEE");
		core::pose::make_pose_from_sequence(pose2, pose2_AND_ss, "fa_standard");

		//Let's make a fauxmetric fauxhexamer!  We need something big enough to have many chains
		core::pose::Pose combo_pose = pose1;
		combo_pose.append_pose_by_jump(pose2, 1);
		combo_pose.append_pose_by_jump(pose1, 1);
		combo_pose.append_pose_by_jump(pose2, 1);
		combo_pose.append_pose_by_jump(pose1, 1);
		combo_pose.append_pose_by_jump(pose2, 1);

		std::string const combo_pose_SS(pose1_AND_ss + pose2_AND_ss + pose1_AND_ss + pose2_AND_ss + pose1_AND_ss + pose2_AND_ss);
		TS_ASSERT_EQUALS(combo_pose_SS.size(), combo_pose.size()); //this is really just assert, but why not make it soft?
		//this is dumb but pose has no string function for this
		for ( core::Size i(1); i<combo_pose.size(); ++i ) {
			combo_pose.set_secstruct(i, combo_pose_SS[i-1]);
		}

		core_init_with_additional_options( "-out:file:output_secondary_structure -out:file:do_not_autoassign_SS");
		core::io::StructFileRepOptionsOP options_op(new core::io::StructFileRepOptions);

		std::ostringstream pdb_stream;
		core::io::pdb::dump_pdb( combo_pose, pdb_stream, options_op);
		//std::cout << pdb_stream.str() << std::endl;

		std::istringstream pdb_result(pdb_stream.str());
		std::string test_answer("");
		//convert the first 66 lines of the PDB into our string to compare to
		for ( core::Size lineNum(1); lineNum <= 66 /*66 SS elements*/; ++lineNum ) {
			std::string templine;
			std::getline(pdb_result, templine);
			test_answer += templine + "\n";
		}

		//commented oddly to keep whitespace intact
		/*HELIX    1   1 HIS A    3  HIS A    8  1                                   6    */
		/*HELIX    2   2 HIS A   11  HIS A   16  1                                   6    */
		/*HELIX    3   3 HIS A   50  HIS A   52  1                                   3    */
		/*HELIX    4   4 HIS A   54  HIS A   55  1                                   2    */
		/*HELIX    5   5 HIS A   57  HIS A   59  1                                   3    */
		/*HELIX    6   6 HIS A   61  HIS A   68  1                                   8    */
		/*HELIX    7   7 HIS A   74  HIS A   75  1                                   2    */
		/*HELIX    8   8 HIS B   80  HIS B   80  1                                   1    */
		/*HELIX    9   9 HIS B  101  HIS B  104  1                                   4    */
		/*HELIX   10  10 HIS B  109  HIS B  112  1                                   4    */
		/*HELIX   11  11 HIS C  139  HIS C  144  1                                   6    */
		/*HELIX   12  12 HIS C  147  HIS C  152  1                                   6    */
		/*HELIX   13  13 HIS C  186  HIS C  188  1                                   3    */
		/*HELIX   14  14 HIS C  190  HIS C  191  1                                   2    */
		/*HELIX   15  15 HIS C  193  HIS C  195  1                                   3    */
		/*HELIX   16  16 HIS C  197  HIS C  204  1                                   8    */
		/*HELIX   17  17 HIS C  210  HIS C  211  1                                   2    */
		/*HELIX   18  18 HIS D  216  HIS D  216  1                                   1    */
		/*HELIX   19  19 HIS D  237  HIS D  240  1                                   4    */
		/*HELIX   20  20 HIS D  245  HIS D  248  1                                   4    */
		/*HELIX   21  21 HIS E  275  HIS E  280  1                                   6    */
		/*HELIX   22  22 HIS E  283  HIS E  288  1                                   6    */
		/*HELIX   23  23 HIS E  322  HIS E  324  1                                   3    */
		/*HELIX   24  24 HIS E  326  HIS E  327  1                                   2    */
		/*HELIX   25  25 HIS E  329  HIS E  331  1                                   3    */
		/*HELIX   26  26 HIS E  333  HIS E  340  1                                   8    */
		/*HELIX   27  27 HIS E  346  HIS E  347  1                                   2    */
		/*HELIX   28  28 HIS F  352  HIS F  352  1                                   1    */
		/*HELIX   29  29 HIS F  373  HIS F  376  1                                   4    */
		/*HELIX   30  30 HIS F  381  HIS F  384  1                                   4    */
		/*SHEET    1   1 1 GLU A  19  GLU A  24  0                                        */
		/*SHEET    1   2 1 GLU A  29  GLU A  35  0                                        */
		/*SHEET    1   3 1 GLU A  37  GLU A  42  0                                        */
		/*SHEET    1   4 1 GLU A  44  GLU A  46  0                                        */
		/*SHEET    1   5 1 GLU A  48  GLU A  48  0                                        */
		/*SHEET    1   6 1 GLU A  71  GLU A  73  0                                        */
		/*SHEET    1   7 1 GLU B  81  GLU B  83  0                                        */
		/*SHEET    1   8 1 GLU B  90  GLU B  95  0                                        */
		/*SHEET    1   9 1 GLU B 116  GLU B 118  0                                        */
		/*SHEET    1  10 1 GLU B 121  GLU B 124  0                                        */
		/*SHEET    1  11 1 GLU B 128  GLU B 131  0                                        */
		/*SHEET    1  12 1 GLU B 135  GLU B 136  0                                        */
		/*SHEET    1  13 1 GLU C 155  GLU C 160  0                                        */
		/*SHEET    1  14 1 GLU C 165  GLU C 171  0                                        */
		/*SHEET    1  15 1 GLU C 173  GLU C 178  0                                        */
		/*SHEET    1  16 1 GLU C 180  GLU C 182  0                                        */
		/*SHEET    1  17 1 GLU C 184  GLU C 184  0                                        */
		/*SHEET    1  18 1 GLU C 207  GLU C 209  0                                        */
		/*SHEET    1  19 1 GLU D 217  GLU D 219  0                                        */
		/*SHEET    1  20 1 GLU D 226  GLU D 231  0                                        */
		/*SHEET    1  21 1 GLU D 252  GLU D 254  0                                        */
		/*SHEET    1  22 1 GLU D 257  GLU D 260  0                                        */
		/*SHEET    1  23 1 GLU D 264  GLU D 267  0                                        */
		/*SHEET    1  24 1 GLU D 271  GLU D 272  0                                        */
		/*SHEET    1  25 1 GLU E 291  GLU E 296  0                                        */
		/*SHEET    1  26 1 GLU E 301  GLU E 307  0                                        */
		/*SHEET    1  27 1 GLU E 309  GLU E 314  0                                        */
		/*SHEET    1  28 1 GLU E 316  GLU E 318  0                                        */
		/*SHEET    1  29 1 GLU E 320  GLU E 320  0                                        */
		/*SHEET    1  30 1 GLU E 343  GLU E 345  0                                        */
		/*SHEET    1  31 1 GLU F 353  GLU F 355  0                                        */
		/*SHEET    1  32 1 GLU F 362  GLU F 367  0                                        */
		/*SHEET    1  33 1 GLU F 388  GLU F 390  0                                        */
		/*SHEET    1  34 1 GLU F 393  GLU F 396  0                                        */
		/*SHEET    1  35 1 GLU F 400  GLU F 403  0                                        */
		/*SHEET    1  36 1 GLU F 407  GLU F 407  0                                        */

		std::string const right_answer = "HELIX    1   1 HIS A    3  HIS A    8  1                                   6    \nHELIX    2   2 HIS A   11  HIS A   16  1                                   6    \nHELIX    3   3 HIS A   50  HIS A   52  1                                   3    \nHELIX    4   4 HIS A   54  HIS A   55  1                                   2    \nHELIX    5   5 HIS A   57  HIS A   59  1                                   3    \nHELIX    6   6 HIS A   61  HIS A   68  1                                   8    \nHELIX    7   7 HIS A   74  HIS A   75  1                                   2    \nHELIX    8   8 HIS B   80  HIS B   80  1                                   1    \nHELIX    9   9 HIS B  101  HIS B  104  1                                   4    \nHELIX   10  10 HIS B  109  HIS B  112  1                                   4    \nHELIX   11  11 HIS C  139  HIS C  144  1                                   6    \nHELIX   12  12 HIS C  147  HIS C  152  1                                   6    \nHELIX   13  13 HIS C  186  HIS C  188  1                                   3    \nHELIX   14  14 HIS C  190  HIS C  191  1                                   2    \nHELIX   15  15 HIS C  193  HIS C  195  1                                   3    \nHELIX   16  16 HIS C  197  HIS C  204  1                                   8    \nHELIX   17  17 HIS C  210  HIS C  211  1                                   2    \nHELIX   18  18 HIS D  216  HIS D  216  1                                   1    \nHELIX   19  19 HIS D  237  HIS D  240  1                                   4    \nHELIX   20  20 HIS D  245  HIS D  248  1                                   4    \nHELIX   21  21 HIS E  275  HIS E  280  1                                   6    \nHELIX   22  22 HIS E  283  HIS E  288  1                                   6    \nHELIX   23  23 HIS E  322  HIS E  324  1                                   3    \nHELIX   24  24 HIS E  326  HIS E  327  1                                   2    \nHELIX   25  25 HIS E  329  HIS E  331  1                                   3    \nHELIX   26  26 HIS E  333  HIS E  340  1                                   8    \nHELIX   27  27 HIS E  346  HIS E  347  1                                   2    \nHELIX   28  28 HIS F  352  HIS F  352  1                                   1    \nHELIX   29  29 HIS F  373  HIS F  376  1                                   4    \nHELIX   30  30 HIS F  381  HIS F  384  1                                   4    \nSHEET    1   1 1 GLU A  19  GLU A  24  0                                        \nSHEET    1   2 1 GLU A  29  GLU A  35  0                                        \nSHEET    1   3 1 GLU A  37  GLU A  42  0                                        \nSHEET    1   4 1 GLU A  44  GLU A  46  0                                        \nSHEET    1   5 1 GLU A  48  GLU A  48  0                                        \nSHEET    1   6 1 GLU A  71  GLU A  73  0                                        \nSHEET    1   7 1 GLU B  81  GLU B  83  0                                        \nSHEET    1   8 1 GLU B  90  GLU B  95  0                                        \nSHEET    1   9 1 GLU B 116  GLU B 118  0                                        \nSHEET    1  10 1 GLU B 121  GLU B 124  0                                        \nSHEET    1  11 1 GLU B 128  GLU B 131  0                                        \nSHEET    1  12 1 GLU B 135  GLU B 136  0                                        \nSHEET    1  13 1 GLU C 155  GLU C 160  0                                        \nSHEET    1  14 1 GLU C 165  GLU C 171  0                                        \nSHEET    1  15 1 GLU C 173  GLU C 178  0                                        \nSHEET    1  16 1 GLU C 180  GLU C 182  0                                        \nSHEET    1  17 1 GLU C 184  GLU C 184  0                                        \nSHEET    1  18 1 GLU C 207  GLU C 209  0                                        \nSHEET    1  19 1 GLU D 217  GLU D 219  0                                        \nSHEET    1  20 1 GLU D 226  GLU D 231  0                                        \nSHEET    1  21 1 GLU D 252  GLU D 254  0                                        \nSHEET    1  22 1 GLU D 257  GLU D 260  0                                        \nSHEET    1  23 1 GLU D 264  GLU D 267  0                                        \nSHEET    1  24 1 GLU D 271  GLU D 272  0                                        \nSHEET    1  25 1 GLU E 291  GLU E 296  0                                        \nSHEET    1  26 1 GLU E 301  GLU E 307  0                                        \nSHEET    1  27 1 GLU E 309  GLU E 314  0                                        \nSHEET    1  28 1 GLU E 316  GLU E 318  0                                        \nSHEET    1  29 1 GLU E 320  GLU E 320  0                                        \nSHEET    1  30 1 GLU E 343  GLU E 345  0                                        \nSHEET    1  31 1 GLU F 353  GLU F 355  0                                        \nSHEET    1  32 1 GLU F 362  GLU F 367  0                                        \nSHEET    1  33 1 GLU F 388  GLU F 390  0                                        \nSHEET    1  34 1 GLU F 393  GLU F 396  0                                        \nSHEET    1  35 1 GLU F 400  GLU F 403  0                                        \nSHEET    1  36 1 GLU F 407  GLU F 407  0                                        \n";

		TS_ASSERT_EQUALS(test_answer, right_answer);

	}

};  // class StructFileRepTests
