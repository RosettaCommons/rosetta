// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/devel/denovo_design/SheetDB.cxxtest.hh
/// @brief  test suite for devel::denovo_design::components::SheetDB
/// @author Tom Linsky (tlinsky@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Unit headers
#include <protocols/denovo_design/components/SheetDB.hh>

// Protocol headers
#include <protocols/denovo_design/architects/StrandArchitect.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/util.hh>

// Project headers
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/fldsgn/topology/util.hh>
#include <protocols/moves/DsspMover.hh>

// Core headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/util/SwitchResidueTypeSet.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Boost headers
#include <boost/assign.hpp>
#include <boost/lexical_cast.hpp>

// C++ headers

#include <devel/denovo_design/components/util.hh>

static THREAD_LOCAL basic::Tracer TR("devel.denovo_design.SheetDB.cxxtest");

// --------------- Test Class --------------- //
class SheetDBTests : public CxxTest::TestSuite {
	// scorefunction
	core::scoring::ScoreFunctionOP scorefxn;

public:

	// Shared initialization goes here.
	void setUp() {
		// load params for ligand
		protocols_init();

		// set preserve header always
		basic::options::option[basic::options::OptionKeys::run::preserve_header].value(true);

		// initialize common filters/movers/scorefxns
		scorefxn = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction() );
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_add_sheet()
	{
		using namespace core::chemical;
		utility::vector1< std::string > params_files;
		ResidueTypeSetCOP const_residue_set = ResidueTypeSetCOP( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
		ResidueTypeSet & residue_set = const_cast< ResidueTypeSet & >( *const_residue_set );
		if ( !residue_set.has_name("D2I") ) params_files.push_back("devel/denovo_design/D2I.params");
		residue_set.read_files_for_custom_residue_types(params_files);

		using namespace protocols::denovo_design::components;

		// create a pose and add it
		core::pose::PoseOP testpose = core::pose::PoseOP( new core::pose::Pose() );
		core::io::pdb::build_pose_from_pdb_as_is( *testpose, "devel/denovo_design/helix1.pdb" );
		TS_ASSERT_EQUALS( testpose->total_residue(), 14 );

		SheetDB db;
		db.set_idealize( false );
		Lengths const lengths = boost::assign::list_of (5)(6)(6);
		StrandOrientations const orientations = boost::assign::list_of (UP)(DOWN)(UP);
		RegisterShifts const register_shifts = boost::assign::list_of (0)(-1)(0);

		core::Size const nstrands = 3;
		std::string const orientation_str = "101";
		std::string const lengths_str = "5:0,6:-1,6:0";
		db.add_sheet( testpose, nstrands, orientation_str, lengths_str, false );

		// map should be saved and retreivable
		SheetList const & poses = db.sheet_list( lengths, orientations, register_shifts );
		TS_ASSERT_EQUALS( poses.size(), 1 );
		TS_ASSERT_EQUALS( poses[1]->total_residue(), testpose->total_residue() );
		TS_ASSERT( poses[1]->is_centroid() );

		// now add a second
		core::pose::PoseOP test2 = core::pose::PoseOP( new core::pose::Pose() );
		core::io::pdb::build_pose_from_pdb_as_is( *test2, "devel/denovo_design/helix15.pdb" );
		TS_ASSERT_EQUALS( test2->total_residue(), 17 );
		db.add_sheet( test2, nstrands, orientation_str, lengths_str, false );

		// now there should be two in the list
		SheetList const & poses2 = db.sheet_list( lengths, orientations, register_shifts );
		TS_ASSERT_EQUALS( poses2.size(), 2 );
		TS_ASSERT_EQUALS( poses2[1]->total_residue(), testpose->total_residue() );
		TS_ASSERT_EQUALS( poses2[2]->total_residue(), test2->total_residue() );
		TS_ASSERT( poses2[1]->is_centroid() );
		TS_ASSERT( poses2[2]->is_centroid() );

		// now add test2 as a different length
		std::string const new_lengths_str = "5:0,5:0,5:0";
		Lengths const new_lengths = boost::assign::list_of(5)(5)(5);
		RegisterShifts const new_shifts = boost::assign::list_of(0)(0)(0);
		db.add_sheet( testpose, nstrands, orientation_str, new_lengths_str, false );
		SheetList const & p3 = db.sheet_list( new_lengths, orientations, new_shifts );
		TS_ASSERT_EQUALS( p3.size(), 1 );
		TS_ASSERT_EQUALS( p3[1]->total_residue(), testpose->total_residue() );
		TS_ASSERT( p3[1]->is_centroid() );
	}

	void test_extract_sheets()
	{
		using namespace protocols::fldsgn::topology;
		using namespace protocols::denovo_design::components;
		utility::vector1< core::pose::PoseOP > poses =
			devel::denovo_design::components::build_poselist_from_multimodel_pdb( "devel/denovo_design/candidate_sheets.multi.pdb.gz" );
		utility::vector1< core::pose::PoseOP > sheets;
		for ( core::Size i=1; i<=poses.size(); ++i ) {
			TR << "Converting pose " << i << std::endl;
			protocols::moves::DsspMover dssp;
			dssp.apply( *poses[i] );
			SS_Info2_OP ssinfo = SS_Info2_OP( new SS_Info2( *poses[i] ) );
			StrandPairingSet spairs = calc_strand_pairing_set( *poses[i], ssinfo );
			TR << spairs << std::endl;
			std::set< core::Size > strands_added;
			core::pose::PoseOP newpose = core::pose::PoseOP( new core::pose::Pose() );
			for ( StrandPairings::const_iterator pair=spairs.begin(); pair != spairs.end(); ++pair ) {
				if ( strands_added.find( (*pair)->s1() ) == strands_added.end() ) {
					add_to_pose( newpose, *poses[i], (*pair)->begin1(), (*pair)->end1() );
					strands_added.insert( (*pair)->s1() );
				}
				if ( strands_added.find( (*pair)->s2() ) == strands_added.end() ) {
					add_to_pose( newpose, *poses[i], (*pair)->begin2(), (*pair)->end2() );
					strands_added.insert( (*pair)->s2() );
				}
			}
		}
	}

	void test_extract_sheet_from_pose()
	{
		using namespace protocols::denovo_design::components;
		core::pose::Pose rsmn;
		core::io::pdb::build_pose_from_pdb_as_is( rsmn, "devel/denovo_design/test_input.pdb" );
		TS_ASSERT_EQUALS( rsmn.size(), 93 );
		protocols::moves::DsspMover dssp;
		dssp.apply( rsmn );

		SheetDB db;
		db.set_idealize( false );
		protocols::fldsgn::topology::SS_Info2_OP ss_info =
			protocols::fldsgn::topology::SS_Info2_OP( new protocols::fldsgn::topology::SS_Info2(rsmn) );
		protocols::fldsgn::topology::StrandPairingSet spairset = protocols::fldsgn::topology::calc_strand_pairing_set( rsmn, ss_info );
		utility::vector1< core::pose::PoseOP > sheets = extract_sheets_from_pose( rsmn, spairset.strand_pairings(), *ss_info, false );

		TS_ASSERT_EQUALS( sheets.size(), 1 );
		core::pose::PoseOP sheet = sheets[1];
		TS_ASSERT_EQUALS( sheet->num_jump(), 2 );
		TS_ASSERT_EQUALS( sheet->conformation().num_chains(), 3 );
		TS_ASSERT_EQUALS( sheet->total_residue(), 6+7+7 );
	}

	void test_extract_sheetset_from_pose()
	{
		using namespace protocols::denovo_design::components;
		core::pose::Pose pose;
		core::io::pdb::build_pose_from_pdb_as_is( pose, "devel/denovo_design/1lvm.pdb" );
		core::util::switch_to_residue_type_set( pose, "centroid" );
		TS_ASSERT_EQUALS( pose.size(), 467 );
		protocols::moves::DsspMover dssp;
		dssp.apply( pose );

		SheetDB db;
		db.set_idealize( false );
		protocols::fldsgn::topology::SS_Info2_OP ss_info =
			protocols::fldsgn::topology::SS_Info2_OP( new protocols::fldsgn::topology::SS_Info2(pose) );
		protocols::fldsgn::topology::StrandPairingSet spairset = protocols::fldsgn::topology::calc_strand_pairing_set( pose, ss_info );
		utility::vector1< core::pose::PoseOP > sheets = extract_sheets_from_pose( pose, spairset.strand_pairings(), *ss_info, false );
		TS_ASSERT_EQUALS( sheets.size(), 8 );
		core::Size max_chains_idx = 0;
		int max_chains = -1;
		for ( core::Size i=1; i<=sheets.size(); ++i ) {
			if ( static_cast<int>(sheets[i]->conformation().num_chains()) > max_chains ) {
				max_chains = static_cast<int>(sheets[i]->conformation().num_chains());
				max_chains_idx = i;
			}
		}
		TR << "max_chains_idx=" << max_chains_idx << std::endl;
		utility::vector1< core::pose::PoseOP > strands = sheets[max_chains_idx]->split_by_chain();
		// should be 7 strands
		TS_ASSERT_EQUALS( static_cast<int>(strands.size()), max_chains );
		TS_ASSERT_EQUALS( strands.size(), 7 );

		for ( core::Size nstrands = 2; nstrands <= strands.size(); ++nstrands ) {
			// for N strands, take all possible N-strand sequential combos
			utility::vector1< core::pose::PoseOP > poses = extract_sheets_from_strandlist( strands, nstrands );
			TR << "nstrands=" << nstrands << " num poses=" << poses.size() << std::endl;
			TS_ASSERT_EQUALS( poses.size(), strands.size() - nstrands + 1 );
			for ( core::Size i=1; i<=poses.size(); ++i ) {
				TS_ASSERT_EQUALS( poses[i]->conformation().num_chains(), nstrands );
				TS_ASSERT_EQUALS( poses[i]->fold_tree().num_jump()+1, nstrands );
			}
		}
	}

	void test_load_sheetlist_from_db()
	{
		using namespace protocols::denovo_design::components;
		core::pose::Pose pose;
		// here we will load a 4-strand, parallel sheet with all lengths=5
		SheetDB db;
		db.set_db_path( "/work/tlinsky/sheet_db/clustered" );

		Lengths const lengths = boost::assign::list_of (5)(5)(5)(5);
		StrandOrientations const orientations = boost::assign::list_of (UP)(UP)(UP)(UP);
		RegisterShifts const shifts = boost::assign::list_of (0)(0)(0)(0);
		SheetList list = db.sheet_list( lengths, orientations, shifts );
		TR << " Loaded " << list.size() << " sheets." << std::endl;
		TS_ASSERT( list.size() );
		SheetList list2 = db.sheet_list_const( lengths, orientations, shifts );
		TS_ASSERT( list2.size() );
		for ( SheetList::iterator s=list.begin(); s!=list.end(); ++s ) {
			TS_ASSERT_EQUALS( (*s)->fold_tree().num_jump()+1, lengths.size() );
			TS_ASSERT_EQUALS( (*s)->conformation().num_chains(), lengths.size() );
			// dssp
			protocols::moves::DsspMover dssp;
			core::pose::PoseOP sh = (*s)->clone();
			dssp.apply( *sh );
			// determine pairings
			protocols::fldsgn::topology::SS_Info2_OP ss_info =
				protocols::fldsgn::topology::SS_Info2_OP( new protocols::fldsgn::topology::SS_Info2(*sh) );
			protocols::fldsgn::topology::StrandPairingSet spairset = protocols::fldsgn::topology::calc_strand_pairing_set( *sh, ss_info );
			for ( core::Size i=1; i<lengths.size(); ++i ) {
				protocols::fldsgn::topology::StrandPairingOP pair = spairset.strand_pairing( i, i+1 );
				TS_ASSERT( pair );
				TS_ASSERT_EQUALS( pair->orient(), 'P' );
				TS_ASSERT_DELTA( pair->rgstr_shift(), 0, 1.0 );
			}

		}
	}

	void test_canonicalize()
	{
		using namespace protocols::denovo_design::components;
		SheetDB db;
		db.set_idealize( false );
		db.set_db_path( "/work/tlinsky/sheet_db/clustered" );
		std::string orientations = "1111";
		std::string lengths = "4:0,6:0,4:2,4:0";
		// Trivial test -- reversing 1111 gives 1111
		std::string rev_orientations = reverse_orientations( orientations );
		std::string rev_lengths = reverse_lengths( orientations, lengths );
		TS_ASSERT_EQUALS( orientations, rev_orientations );
		TS_ASSERT_EQUALS( rev_lengths, "4:0,4:0,6:-2,4:0" );
		// double reverse should give back original
		TS_ASSERT_EQUALS( orientations, reverse_orientations( rev_orientations ) );
		TS_ASSERT_EQUALS( lengths, reverse_lengths( rev_orientations, rev_lengths ) );
		TS_ASSERT_EQUALS( choose_canonical( lengths, rev_lengths ), choose_canonical( rev_lengths, lengths ) );
		TS_ASSERT_EQUALS( canonicalize( orientations, lengths ), canonicalize( rev_orientations, rev_lengths ) );

		// now try antiparallel
		orientations = "1110";
		rev_orientations = reverse_orientations( orientations );
		rev_lengths = reverse_lengths( orientations, lengths );
		// in this case, orientations are identical
		TS_ASSERT_EQUALS( rev_orientations, "1000" );
		TS_ASSERT_EQUALS( rev_lengths, "4:0,4:0,6:0,4:2" );

		// double reverse should again give back original
		TS_ASSERT_EQUALS( orientations, reverse_orientations( rev_orientations ) );
		TS_ASSERT_EQUALS( lengths, reverse_lengths( rev_orientations, rev_lengths ) );
		TS_ASSERT_EQUALS( choose_canonical( lengths, rev_lengths ), choose_canonical( rev_lengths, lengths ) );
		TS_ASSERT_EQUALS( canonicalize( orientations, lengths), canonicalize( rev_orientations, rev_lengths) );
	}


	void test_wrongorder_in_db() {
		using namespace protocols::denovo_design::components;

		SheetDB sheetdb;
		sheetdb.set_db_path( "/work/tlinsky/sheet_db/clustered" );
		core::Size const nstrands = 3;
		std::string const orientations_str = "101";
		std::string const lengths_str = "7:0,7:0,5:0";
		Lengths const lengths = boost::assign::list_of (7)(7)(5);
		StrandOrientations const orientations = boost::assign::list_of (UP)(DOWN)(UP);
		RegisterShifts const shifts = boost::assign::list_of (0)(0)(0);

		std::string const id = "testsheet";
		std::string const ids[] = { "s1", "s2", "s3" };
		core::Size const len[] = { 7, 7, 5 };

		// construct permutation
		StructureDataOP perm( new StructureData( "testsheet" ) );
		for ( core::Size i=1; i<=nstrands; ++i ) {
			std::string ss = "L";
			utility::vector1< std::string > abego;
			abego.push_back( "X" );
			for ( core::Size r=1; r<=len[i-1]; ++r ) {
				ss += "E";
				abego.push_back( "B" );
			}
			ss += "L";
			abego.push_back( "X" );
			StructureData subperm( "test1" );
			subperm.add_segment( ids[i-1],
				Segment( ss, protocols::denovo_design::abego_str(abego), false, false ) );
			perm->merge( subperm );
		}
		perm->add_prefix_to_segments( id );

		core::Size count = 0;
		core::Size failcount = 0;
		SheetList const & list = sheetdb.sheet_list( lengths, orientations, shifts );
		for ( SheetList::const_iterator s=list.begin(), end=list.end(); s != end; ++s ) {
			++count;
			core::pose::PoseCOP chosensheet = *s;
			bool failed = false;
			// check for improper terminal variant
			for ( int i=1, end=nstrands; i<=end; ++i ) {
				std::string const strandid = ids[i-1];
				Segment const & r = perm->segment( id + "." + strandid );
				for ( core::Size resnum=r.lower()+1; resnum<r.upper(); ++resnum ) {
					TS_ASSERT( !chosensheet->residue(resnum).is_terminus() );
					if ( chosensheet->residue(resnum).is_terminus() ) {
						failed = true;
						TR << "Residue " << resnum << " has an incorrect terminal variant in the middle of a component."  << std::endl;
						TR << "StructureData=" << *perm << std::endl;
						for ( core::Size j=1; j<=chosensheet->total_residue(); ++j ) {
							TR << "RES " << chosensheet->residue(j).name() << j << " : " << chosensheet->chain(j) << std::endl;
						}
						chosensheet->dump_pdb( "bad_terminus.pdb" );
					}
				}
				if ( failed ) {
					TR << "number " << count << " failed." << std::endl;
					break;
				}
			}
			if ( failed ) {
				++failcount;
			}
		}
		TR << failcount << " of " << count << " failed." << std::endl;
	}

};
