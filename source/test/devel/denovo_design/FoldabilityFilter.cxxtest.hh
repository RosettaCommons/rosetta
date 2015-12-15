// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/devel/denovo_design/FoldabilityFilter.cxxtest.hh
/// @brief  test suite for devel::denovo_design::components::FoldabilityFilter
/// @author Tom Linsky (tlinsky@uw.edu)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Unit headers
#include <devel/denovo_design/filters/FoldabilityFilter.hh>

// Project headers

// Protocol headers
#include <protocols/matdes/SymDofMover.hh>
#include <protocols/moves/DsspMover.hh>

// Core headers
#include <core/kinematics/FoldTree.hh>
#include <core/io/pdb/file_data.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/ScoreFunction.hh>

// Utility headers
#include <basic/datacache/DataMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// C++ headers
#include <test/protocols/denovo_design/test_utils.hh>
static THREAD_LOCAL basic::Tracer TR("devel.denovo_design.FoldabilityFilter.cxxtest");

// Dummy class so that I can access protected functions
class ProtectedFoldabilityFilter : public devel::denovo_design::filters::FoldabilityFilter {
public:
	/// @brief gets aa string, ss string, and abego vector for the area to rebuild
	void Tget_aa_ss_abego(
		std::string & aa,
		std::string & ss,
		utility::vector1< std::string > & abego,
		core::Size const start,
		core::Size & end,
		core::pose::Pose const & pose ) const
	{
		get_aa_ss_abego( aa, ss, abego, start, end, pose );
	}

	core::pose::PoseOP Tgenerate_pose( core::pose::Pose const & pose ) const
	{
		return generate_pose( pose );
	}

	/// @brief deletes the segment from start to end (inclusive) from the pose
	void Tprepare_pose(
		core::pose::Pose & pose,
		core::Size const start,
		core::Size const end ) const
	{
		prepare_pose( pose, start, end );
	}

	/// @brief abegod scores
	core::Real Tabegodb_score(
		core::pose::Pose const & pose,
		utility::vector1< std::string > const & abego,
		core::Size const start,
		core::Size const stop ) const
	{
		return abegodb_score( pose, abego, start, stop );
	}

};

// --------------- Test Class --------------- //
class FoldabilityFilterTests : public CxxTest::TestSuite {
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

	void test_foldability() {
		using namespace devel::denovo_design::filters;
		core::pose::Pose input_pose;
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "devel/denovo_design/test_foldability.pdb" );
		// Add symmetry
		std::stringstream stream;
		stream << "<SymDofMover name=gen_docked_config symm_file=\"devel/denovo_design/C5_Z.sym\" sym_dof_names=\"JS1\" />" << std::endl;
		utility::tag::TagCOP symtag = utility::tag::Tag::create( stream );
		basic::datacache::DataMap data;
		protocols::filters::Filters_map filters;
		protocols::moves::Movers_map movers;
		protocols::matdes::SymDofMover sym;
		sym.parse_my_tag( symtag, data, filters, movers, input_pose );
		sym.apply( input_pose );
		TS_ASSERT( core::pose::symmetry::is_symmetric(input_pose) );

		ProtectedFoldabilityFilter folder;
		folder.add_segment( 5, 15 );

		// test pose copy generation -- should no longer be symmetric
		core::pose::PoseOP posecopy = folder.Tgenerate_pose( input_pose );
		TS_ASSERT( !core::pose::symmetry::is_symmetric(*posecopy) );
		TS_ASSERT_EQUALS( posecopy->total_residue(), 371 );
		for ( core::Size i=1, endi=posecopy->total_residue(); i<=endi; ++i ) {
			TS_ASSERT( posecopy->residue_type(i).name() != "VRT" );
			TR << i << " : " << posecopy->residue(i).name() << std::endl;
		}

		// test generating start/end residues
		core::Size start = 5;
		core::Size end = 15;

		// dssp the pose
		protocols::moves::DsspMover dssp;
		dssp.apply( *posecopy );

		// test generating aa/ss/abego strings
		std::string aa = "";
		std::string ss = "";
		utility::vector1< std::string > abego;
		folder.Tget_aa_ss_abego( aa, ss, abego, start, end, *posecopy );
		TS_ASSERT_EQUALS( aa.size(), posecopy->total_residue()-1 );
		TS_ASSERT_EQUALS( ss.size(), posecopy->total_residue()-1 );
		TS_ASSERT_EQUALS( abego.size(), posecopy->total_residue()-1 );
		TR << "AA=" << aa << std::endl;
		TR << "SS=" << ss << std::endl;

		// test deleting the segment in question
		folder.Tprepare_pose( *posecopy, start, end );
		// should be smaller by 1 for the end+1 residue that should be cut
		TS_ASSERT_EQUALS( posecopy->total_residue(), 370 );
		// should be one jump
		TS_ASSERT_EQUALS( posecopy->fold_tree().num_jump(), 1 );
		// that jump should start at residue start-1
		// it should end at residue start since the residues in the segment have been deleted
		TS_ASSERT_EQUALS( posecopy->fold_tree().jump_edge(1).start(), start - 1 );
		TS_ASSERT_EQUALS( posecopy->fold_tree().jump_edge(1).stop(), end + 1 );
	}

};
