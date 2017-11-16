// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/devel/denovo_design/AlignResiduesMoverTests.cxxtest.hh
/// @brief  test suite for devel::denovo_design::components::AlignResiduesMover
/// @author Tom Linsky (tlinsky@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Unit headers
#include <protocols/denovo_design/movers/AlignResiduesMover.hh>
#include <protocols/denovo_design/components/ExtendedPoseBuilder.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/residue_selectors/NamedSegmentSelector.hh>

// Protocol headers

// Project headers

// Core headers
#include <core/conformation/Residue.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <test/protocols/denovo_design/test_utils.hh>
static basic::Tracer TR("devel.denovo_design.AlignResiduesMoverTests.cxxtest");

using namespace protocols::denovo_design::components;
using namespace protocols::denovo_design::movers;

// --------------- Test Class --------------- //
class AlignResiduesMoverTests : public CxxTest::TestSuite {
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

	void test_tomponent_cstfile()
	{
		/*
		TomponentFileConstraintGenerator myrcg;
		myrcg.set_cstfile( "devel/denovo_design/test.cstfile" );
		StructureDataOP perm( new StructureData("unittest") );
		utility::vector1< std::string > abegos;
		abegos.push_back("X");
		for ( int i=1; i<=8; ++i ) {
		abegos.push_back("A");
		}
		abegos.push_back("X");

		std::string const name1 = "build_dyad.ser_his.1";
		std::string const name2 = "build_dyad.ser_his.2";
		StructureData subperm1( name1 );
		subperm1.add_segment( name1, Segment( "LHHHHHHHHL", abegos, 1, false, false ) );
		StructureData subperm2( name2 );
		subperm2.add_segment( name2, Segment( "LHHHHHHHHL", abegos, 1, false, false ) );
		perm->merge(subperm1);
		perm->merge(subperm2);

		AddAlias add_alias;
		add_alias.set_id( "build_dyad.cat_res.1" );
		add_alias.set_segment( "build_dyad.ser_his.2" );
		add_alias.set_resid( 4 );
		add_alias.setup_and_apply_permutation( *perm );
		TS_ASSERT_EQUALS( add_alias.get_last_move_status(), protocols::moves::MS_SUCCESS );
		TS_ASSERT_EQUALS( perm->alias( "build_dyad.cat_res.1" ), 15 );

		AddAlias add_absolute_alias;
		add_absolute_alias.set_id( "myalias" );
		add_absolute_alias.set_resid( 16 );
		add_absolute_alias.setup_and_apply_permutation( *perm );
		TS_ASSERT_EQUALS( add_absolute_alias.get_last_move_status(), protocols::moves::MS_SUCCESS );
		TS_ASSERT_EQUALS( perm->alias( "myalias" ), 16 );

		std::string const csts_str = myrcg.get_constraint_string( *perm );
		TR << "New csts_str = " << csts_str << std::endl;
		TS_ASSERT_EQUALS( csts_str.find("%%"), std::string::npos );
		utility::vector1< std::string > lines = utility::string_split( csts_str, '\n' );
		for ( core::Size i=1; i<=lines.size(); ++i ) {
		utility::vector1< std::string > fields = utility::string_split( lines[i] );
		if ( i==1 ) {
		TS_ASSERT_EQUALS( fields[1], "AtomPair" );
		TS_ASSERT_EQUALS( fields[2], "O2" );
		TS_ASSERT_EQUALS( fields[3], "2" ); // substitution 1
		TS_ASSERT_EQUALS( fields[4], "N" );
		TS_ASSERT_EQUALS( fields[5], "3" );
		TS_ASSERT_EQUALS( fields[6], "HARMONIC" );
		TS_ASSERT_EQUALS( fields[7], "2.8" );
		TS_ASSERT_EQUALS( fields[8], "0.3" );
		} else if ( i==2 ) {
		TS_ASSERT_EQUALS( fields[1], "Angle" );
		TS_ASSERT_EQUALS( fields[2], "O2" );
		TS_ASSERT_EQUALS( fields[3], "4" ); // substitution 2
		TS_ASSERT_EQUALS( fields[4], "N" );
		TS_ASSERT_EQUALS( fields[5], "15" );
		TS_ASSERT_EQUALS( fields[6], "H" );
		TS_ASSERT_EQUALS( fields[7], "15" );
		TS_ASSERT_EQUALS( fields[8], "HARMONIC" );
		TS_ASSERT_EQUALS( fields[9], "0.0" );
		TS_ASSERT_EQUALS( fields[10], "0.1" );
		} else {
		TS_ASSERT_EQUALS( lines[i].size(), 0 );
		}
		}
		*/
	}

	void test_align_theozyme()
	{
		using namespace protocols::denovo_design::residue_selectors;
		core::pose::Pose pose;

		core::io::pdb::build_pose_from_pdb_as_is( pose, "devel/denovo_design/cat_residues.pdb" );
		core::pose::Pose const input_pose = pose;

		core::pose::Pose ph1, ph2;
		core::io::pdb::build_pose_from_pdb_as_is( ph1, "protocols/denovo_design/components/helix14.pdb" );
		core::io::pdb::build_pose_from_pdb_as_is( ph2, "protocols/denovo_design/components/helix15.pdb" );

		StructureData sd( "prot" );

		Segment cat1( "prot.cat.1", "L", "X", true, true );
		cat1.set_template_pose( input_pose, 3, 3 );
		Segment cat2( "prot.cat.2", "L", "X", true, true );
		cat2.set_template_pose( input_pose, 7, 7 );

		Segment h1( "prot.h1", "LHHHHHHHHHHHHHHL", "XAAAAAAAAAAAAAAX", false, false );
		h1.set_template_pose( ph1, 2, 15 );
		h1.set_movable_group( 2 );

		Segment h2( "prot.h2", "LHHHHHHHHHHHHHHHL", "XAAAAAAAAAAAAAAAX", false, false );
		h2.set_movable_group( 3 );
		h2.set_template_pose( ph2, 2, 16 );

		sd.add_segment( cat1 );
		sd.add_segment( cat2 );
		sd.add_segment( h1 );
		sd.add_segment( h2 );

		ExtendedPoseBuilder builder;
		core::pose::PoseOP pose_ptr = builder.apply( sd );
		TS_ASSERT( pose_ptr );
		core::pose::Pose & built = *pose_ptr;
		TS_ASSERT_THROWS_NOTHING( sd.check_pose_consistency( built ) );

		StructureData const orig = sd;
		core::pose::Pose const orig_pose = built;

		// create alignresidues object
		std::string const myid = "align";
		AlignResiduesMover align;
		align.set_id( myid );

		NamedSegmentSelector selcat1( "prot.cat.1", "1" );
		NamedSegmentSelector selcat2( "prot.cat.2", "1" );
		align.add_template_selector( selcat1 );
		align.add_template_selector( selcat2 );

		NamedSegmentSelector sel( "prot.h1", "4" );
		NamedSegmentSelector sel2( "prot.h2", "8" );
		align.add_target_selector( sel );
		align.add_target_selector( sel2 );

		align.apply( built );
		sd = StructureDataFactory::get_instance()->get_from_const_pose( built );
		sd.check_pose_consistency( built );
		TS_ASSERT_EQUALS( align.get_last_move_status(), protocols::moves::MS_SUCCESS );

		for ( core::Size r=1; r<=built.size(); ++r ) {
			TR << built.residue( r ).name() <<  " " << r << std::endl;
		}

		// aliases should be set
		std::string const alias1 = myid + ".1";
		std::string const alias2 = myid + ".2";

		TS_ASSERT_EQUALS( sd.alias( alias1 ), 5 );
		TS_ASSERT_EQUALS( sd.alias( alias2 ), 25 );

		// h1,h2 in same MG
		TS_ASSERT_EQUALS(
			sd.segment("prot.h1").movable_group(),
			sd.segment("prot.h2").movable_group() );

		// cat.1 should be ser, cat.2 should be his
		TS_ASSERT_EQUALS( built.residue( sd.alias( alias1 ) ).name1(), 'S' );
		TS_ASSERT_EQUALS( built.residue( sd.alias( alias2 ) ).name1(), 'H' );

		// h1 should be aligned with prot.cat.1
		TS_ASSERT_DELTA(
			built.residue( sd.alias( alias1 ) ).xyz("CA").distance(
			orig_pose.residue(orig.segment("prot.cat.1").start()).xyz("CA") ), 0.0, 1e-6 );

		// h2 should be aligned with prot.cat.2
		TS_ASSERT_DELTA(
			built.residue( sd.alias( alias2 ) ).xyz("CA").distance(
			orig_pose.residue(orig.segment("prot.cat.2").start()).xyz("CA") ), 0.0, 1e-6 );

		// cat1 --> cat2 should be same as h1 --> h2
		core::Vector const catvec1 = orig_pose.residue(orig.segment("prot.cat.1").start()).xyz("CA");
		core::Vector const catvec2 = orig_pose.residue(orig.segment("prot.cat.2").start()).xyz("CA");
		core::Vector const catvecs = catvec1 - catvec2;

		core::Vector const hvec1 = built.residue( sd.alias( alias1 ) ).xyz("CA");
		core::Vector const hvec2 = built.residue( sd.alias( alias2 ) ).xyz("CA");
		core::Vector const hvecs = hvec1 - hvec2;

		core::Vector const bcatvec1 = orig_pose.residue(orig.segment("prot.cat.1").start()).xyz("CB");
		core::Vector const bcatvec2 = orig_pose.residue(orig.segment("prot.cat.2").start()).xyz("CB");
		core::Vector const bcatvecs = bcatvec1 - bcatvec2;

		core::Vector const bhvec1 = built.residue( sd.alias( alias1 ) ).xyz("CB");
		core::Vector const bhvec2 = built.residue( sd.alias( alias2 ) ).xyz("CB");
		core::Vector const bhvecs = bhvec1 - bhvec2;

		TS_ASSERT_DELTA( catvecs.length(), hvecs.length(), 1e-6 );
	}

};
