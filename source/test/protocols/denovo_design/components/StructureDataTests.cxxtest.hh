// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/denovo_design/StructureData.cxxtest.hh
/// @brief  test suite for protocols::denovo_design::components::StructureData
/// @author Tom Linsky (tlinsky@uw.edu)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Unit headers
#include <protocols/denovo_design/components/FoldGraph.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>

// Protocol headers
#include <protocols/denovo_design/util.hh>

// Core headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
#include <utility/tag/Tag.hh>

// Boost headers
#include <boost/assign.hpp>

// C++ headers
static THREAD_LOCAL basic::Tracer TR("protocols.denovo_design.StructureData.cxxtest");

// --------------- Test Class --------------- //
class StructureDataTests : public CxxTest::TestSuite {
	typedef protocols::denovo_design::SegmentNameSet SegmentNameSet;
	typedef protocols::denovo_design::SegmentNameList SegmentNameList;
	typedef protocols::denovo_design::components::StructureData StructureData;
	typedef protocols::denovo_design::components::StructureDataOP StructureDataOP;
	typedef protocols::denovo_design::components::StructureDataFactory StructureDataFactory;
public:

	// Shared initialization goes here.
	void setUp() {
		// load params for ligand
		protocols_init_with_additional_options( "-extra_res_fa protocols/denovo_design/QTS.params -extra_patch_fa protocols/denovo_design/QTS_connectC6.txt -extra_patch_fa protocols/denovo_design/CYS_connectSG.txt" );

		// set preserve header always
		basic::options::option[basic::options::OptionKeys::run::preserve_header].value(true);
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_delete_segment()
	{
		using namespace protocols::denovo_design;
		using namespace protocols::denovo_design::components;
		std::ifstream in_xml( "protocols/denovo_design/test_sd.xml" );

		StructureData perm = StructureDataFactory::get_instance()->create_from_xml( in_xml );
		StructureData const orig = perm;

		for ( SegmentNameList::const_iterator c=perm.segments_begin(); c!=perm.segments_end(); ++c ) {
			TS_ASSERT( perm.segment( *c ).stop() >= perm.segment( *c ).start() );
			TS_ASSERT( perm.segment( *c ).safe() >= perm.segment( *c ).start() );
		}

		TR << "ORIG: " << orig << std::endl;
		// delete bb.h1_sheet -- it shouldn't mess up any residues
		perm.delete_segment( "bb.h1_sheet" );
		TR << "DELE: " << perm << std::endl;
		for ( SegmentNameList::const_iterator c=perm.segments_begin(); c!=perm.segments_end(); ++c ) {
			TR.Debug << "Looking at segment " << *c << std::endl;
			TS_ASSERT( orig.has_segment( *c ) );
			TS_ASSERT( perm.segment( *c ).stop() >= perm.segment( *c ).start() );
			TS_ASSERT( perm.segment( *c ).safe() >= perm.segment( *c ).start() );
			TS_ASSERT_EQUALS( perm.segment( *c ).elem_length(), orig.segment( *c ).elem_length() );
			TS_ASSERT_EQUALS( perm.segment( *c ).safe() - perm.segment( *c ).start(),
				orig.segment( *c ).safe() - orig.segment( *c ).start() );
		}
	}

	/// @brief helper function for test_move_segments
	void move_pair( StructureData const & orig, std::string const & seg1, std::string const & seg2 )
	{
		StructureData sd = orig;
		sd.move_segment( seg1, seg2 );

		std::stringstream orig_out;
		orig_out << orig;
		std::stringstream new_out;
		new_out << sd;
		TS_ASSERT_EQUALS( new_out.str(), orig_out.str() );
	}

	void test_move_segments()
	{
		// Create a test StructureData
		std::string const sd_id = "TestCreateLoops";
		std::stringstream xml;
		xml << "<StructureData name=" << sd_id << " length=34 pose_length=38 >" << std::endl;
		xml << "<ResidueRange name=h1 start=2 group=1 nterm=0 cterm=1 ss=LHHHHHHHHHH abego=XAAAAAAAAAA upper_segment=l1 />" << std::endl;
		xml << "<ResidueRange name=l1 start=12 group=3 nterm=1 cterm=1 ss=LLLL abego=BGAB lower_segment=h1 upper_segment=h2 />" << std::endl;
		xml << "<ResidueRange name=h2 start=16 group=2 nterm=1 cterm=0 ss=HHHHHHHHL abego=AAAAAAAAX lower_segment=l1 />" << std::endl;
		xml << "<ResidueRange name=fixed start=25 nterm=0 cterm=0 group=1 ss=LEEEEELLEEEEEL abego=XBBBBBGGBBBBBX />" << std::endl;
		xml << "</StructureData>" << std::endl;
		StructureData sd = StructureDataFactory::get_instance()->create_from_xml( xml );
		StructureData const orig = sd;
		std::stringstream orig_ss;
		orig_ss << orig;

		// move h1 to after fixed -- h1, l1, h2 should all move, and fixed should be first
		sd.move_segment( "fixed", "h1" );
		TS_ASSERT_EQUALS( sd.segment("fixed").lower(), 1 );
		TS_ASSERT_EQUALS( sd.segment("fixed").upper() + 1, sd.segment("h1").lower() );
		TS_ASSERT_EQUALS( sd.segment("h1").upper() + 1, sd.segment("l1").lower() );
		TS_ASSERT_EQUALS( sd.segment("l1").upper() + 1, sd.segment("h2").lower() );

		// check correct connection status
		TS_ASSERT_EQUALS( sd.segment("fixed").cterm_included(), false );
		TS_ASSERT_EQUALS( sd.segment("h1").nterm_included(), false );
		TS_ASSERT_EQUALS( sd.segment("fixed").upper_segment(), "" );
		TS_ASSERT_EQUALS( sd.segment("h1").lower_segment(), "" );

		sd = orig;
		sd.move_segment( "h1", "l1" );

		move_pair( orig, "h1", "l1" );
		move_pair( orig, "l1", "h2" );
		move_pair( orig, "h2", "fixed" );

		// bad combination of segments should throw error
		TS_ASSERT_THROWS_ANYTHING( move_pair( orig, "fixed", "l1" ) );
		TS_ASSERT_THROWS_ANYTHING( move_pair( orig, "fixed", "h2" ) );
		TS_ASSERT_THROWS_ANYTHING( move_pair( orig, "h1", "h2" ) );
		TS_ASSERT_THROWS_ANYTHING( move_pair( orig, "h1", "fixed" ) );
		TS_ASSERT_THROWS_ANYTHING( move_pair( orig, "l1", "h1" ) );
		TS_ASSERT_THROWS_ANYTHING( move_pair( orig, "l1", "fixed" ) );
		TS_ASSERT_THROWS_ANYTHING( move_pair( orig, "h2", "h1" ) );
		TS_ASSERT_THROWS_ANYTHING( move_pair( orig, "h2", "l1" ) );
		TS_ASSERT_THROWS_ANYTHING( move_pair( orig, "h1", "h1" ) );
		TS_ASSERT_THROWS_ANYTHING( move_pair( orig, "l1", "l1" ) );
		TS_ASSERT_THROWS_ANYTHING( move_pair( orig, "h2", "h2" ) );
		TS_ASSERT_THROWS_ANYTHING( move_pair( orig, "fixed", "fixed" ) );
	}

	void test_non_peptidic_bonds()
	{
		using namespace protocols::denovo_design;
		using namespace protocols::denovo_design::components;
		core::pose::Pose input_pose;
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/denovo_design/test_structuredata_nonpeptidic.pdb" );
		StructureData const sd = StructureDataFactory::get_instance()->get_from_pose( input_pose );

		core::conformation::Residue const prev_cys = input_pose.residue( sd.alias( "bb.cys" ) );
		core::conformation::Residue const prev_qts = input_pose.residue( sd.alias( "bb.ligand" ) );

		// Foldgraph should be able to handle non-peptide bonds
		FoldGraph fg( sd );
		core::kinematics::FoldTree const ft = fg.fold_tree( boost::assign::list_of("1") );
		TS_ASSERT( ft.check_fold_tree() );
		input_pose.fold_tree( ft );

		// rotate a chi angle in QTS
		core::Size const cysres = sd.alias( "bb.cys" );
		core::Size const ligres = sd.alias( "bb.ligand" );
		for ( core::Size chi=1; chi<=6; ++chi ) {
			input_pose.set_chi( chi, ligres, 50.0 );
		}

		TS_ASSERT_DELTA( prev_cys.xyz( "CB" ).distance( input_pose.residue( cysres ).xyz( "CB" ) ), 0.0, 1e-4 );
		TS_ASSERT_DELTA( prev_cys.xyz( "SG" ).distance( input_pose.residue( cysres ).xyz( "SG" ) ), 0.0, 1e-4 );
		TS_ASSERT_DELTA( prev_qts.xyz( "C6" ).distance( input_pose.residue( ligres ).xyz( "C6" ) ), 0.0, 1e-4 );
		TS_ASSERT_DELTA( prev_qts.xyz( "O2" ).distance( input_pose.residue( ligres ).xyz( "O2" ) ), 0.0, 1e-4 );

		// rotate chi in cys
		input_pose.set_chi( 1, cysres, 50.0 );

		TS_ASSERT_DELTA( prev_cys.xyz( "CB" ).distance( input_pose.residue( cysres ).xyz( "CB" ) ), 0.0, 1e-4 );
		TS_ASSERT( prev_cys.xyz( "SG" ).distance( input_pose.residue( cysres ).xyz( "SG" ) ) > 1e-4 );
		TS_ASSERT( prev_qts.xyz( "C6" ).distance( input_pose.residue( ligres ).xyz( "C6" ) ) > 1e-4 );
		TS_ASSERT( prev_qts.xyz( "O2" ).distance( input_pose.residue( ligres ).xyz( "O2" ) ) > 1e-4 );
	}

	void test_enzdes_remarks()
	{
		using namespace core::chemical;
		using namespace protocols::denovo_design::components;

		// Residue definitions can't be supplied on the command line b/c
		// the ResidueTypeSet is already initialized.
		utility::vector1< std::string > params_files;
		ResidueTypeSetCOP const_residue_set = ResidueTypeSetCOP( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
		ResidueTypeSet & residue_set = const_cast< ResidueTypeSet & >(*const_residue_set);
		if ( !residue_set.has_name("D2I") ) params_files.push_back("devel/denovo_design/D2I.params");
		residue_set.read_files_for_custom_residue_types(params_files);

		core::pose::Pose input_pose;
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/denovo_design/components/test_enz_remark_input.pdb" );
		TS_ASSERT( input_pose.pdb_info() );

		// Get original remarks
		core::io::Remarks const origremarks = input_pose.pdb_info()->remarks();

		// Create StructureData to fit this pose
		std::string const dummy_ss( 92, 'L' );
		std::string const dummy_abego( 92, 'X' );
		std::stringstream sd_xml;
		sd_xml << "<StructureData name=\"test_enzdes_remarks\" length=\"91\" pose_length=\"93\" >" << std::endl;
		sd_xml << "<ResidueRange name=\"1\" start=\"2\" nterm=\"0\" cterm=\"0\" group=\"1\" ss=\"" << dummy_ss
			<< "\" abego=\"" << dummy_abego << "\" />" << std::endl;
		sd_xml << "<ResidueRange name=\"2\" start=\"93\" nterm=\"1\" cterm=\"1\" group=\"1\" ss=\"L\" abego=\"X\" />" << std::endl;
		sd_xml << "</StructureData>" << std::endl;

		StructureData sd = StructureDataFactory::get_instance()->create_from_xml( sd_xml );

		// to start, there should be no remarks
		TS_ASSERT( sd.retrieve_remarks( input_pose ).empty() );

		// Save remarks into a structure data
		sd.save_remarks( origremarks );
		TR << sd << std::endl;

		// We should then be able to retrieve identical remarks
		core::io::Remarks const newremarks = sd.retrieve_remarks( input_pose );

		TS_ASSERT_EQUALS( origremarks.size(), newremarks.size() );
		for ( core::io::Remarks::const_iterator orig=origremarks.begin(), n=newremarks.begin();
				( orig != origremarks.end() ) && ( n != newremarks.end() ); ++orig, ++n ) {
			TR << "Testing " << *n << " vs " << *orig << std::endl;
			TS_ASSERT_EQUALS( n->num, orig->num );
			TS_ASSERT_EQUALS( n->value, orig->value );
		}
	}

	void test_slice() {
		utility::io::izstream xml( "protocols/denovo_design/test_sd.xml" );
		StructureData sd = StructureDataFactory::get_instance()->create_from_xml( xml );

		SegmentNameSet const slice_segments = boost::assign::list_of("bb.sheet.s1")("bb.sheet_h1")("bb.h1")("bb.h1_sheet")("bb.sheet.s2");
		StructureData const sliced = sd.slice( slice_segments, false );

		for ( SegmentNameList::const_iterator s=sd.segments_begin(); s!=sd.segments_end(); ++s ) {
			if ( slice_segments.find( *s ) == slice_segments.end() ) {
				TS_ASSERT( !sliced.has_segment( *s ) );
			} else {
				TS_ASSERT( sliced.has_segment( *s ) );
			}
		}
	}
};
