// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/calc_taskop_movers/ForceDisulfidesMover.cxxtest.hh
/// @brief  test for ForceDisulfidesMover
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/rosettascripts.hh>
#include <test/util/schema_utilities.hh>
#include <test/util/pose_funcs.hh>

// Project Headers

#include <protocols/calc_taskop_movers/ForceDisulfidesMover.hh>

#include <core/conformation/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/types.hh>

// Utility Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.calc_taskop_movers/ForceDisulfidesMover.cxxtest.hh");

// --------------- Test Class --------------- //

class ForceDisulfidesMoverTests : public CxxTest::TestSuite {

private:
	core::pose::PoseCOP pose_;
	core::scoring::ScoreFunctionOP scorefxn_;
public:

	void setUp() {
		core_init();
		pose_ = core::import_pose::pose_from_file("core/scoring/4dO8B.pdb");
		scorefxn_ = utility::pointer::make_shared< core::scoring::ScoreFunction >();
	}

	void tearDown() {
	}

	void test_parse() {
		using namespace protocols::calc_taskop_movers;
		check_if_mover_tag_validates< ForceDisulfidesMover >("<ForceDisulfides scorefxn=\"beta_nov16\" name=\"disulf\" disulfides=\"8:463,470:474,278:302,46:274,58:70,91:136\" />");
		check_if_mover_tag_validates< ForceDisulfidesMover >("<ForceDisulfides name=\"disulf\" disulfides=\"8A:463A,470B:474B,278C:302C,46B:274C,58H:70L,91B:136A\" />");
		check_if_mover_tag_validates< ForceDisulfidesMover >("<ForceDisulfides name=\"disulf\" disulfides=\"8A:463A,470:474,278C:302B,46:274,58:70H,91L:136\" />");
		check_if_mover_tag_validates< ForceDisulfidesMover >("<ForceDisulfides name=\"disulf\" disulfides=\"8A:463A,470B:474B,278C:302C,46B:274C,58H:70L,91B:136A\" remove_existing=\"true\" />");
		check_if_mover_tag_validates< ForceDisulfidesMover >("<ForceDisulfides name=\"disulf\" remove_existing=\"true\" />");
		check_if_mover_tag_validates< ForceDisulfidesMover >("<ForceDisulfides name=\"disulf\" remove_existing=\"true\" repack=\"false\" />");
	}

	void test_original() {
		// This is mainly just to check that the original pose is constructed as we expect.
		// If this is failing, it's an indication that the input for this unit test is messed up.
		std::string original_sequence = pose_->annotated_sequence();
		TS_ASSERT_EQUALS(original_sequence, "L[LEU:NtermProteinFull]TC[CYS:disulfide]VTSKSIFGITTENC[CYS:disulfide]PDGQNLC[CYS:disulfide]FKKWYYIVPRYSDITWGC[CYS:disulfide]AATC[CYS:disulfide]PKPTNVRETIRC[CYS:disulfide]C[CYS:disulfide]ETDKC[CYS:disulfide]NE[GLU:CtermProteinFull]");

		utility::vector1< std::pair<core::Size,core::Size> > disulfides;
		core::conformation::disulfide_bonds( pose_->conformation(), disulfides);
		TS_ASSERT_EQUALS( disulfides.size(), 4 );
		// We make the somewhat reasonable assumption that the pairs are ordered internally
		// (But we don't assume that any disulfide is in a particular order.
		TS_ASSERT( disulfides.contains( std::make_pair(3,24) ) );
		TS_ASSERT( disulfides.contains( std::make_pair(17,42) ) );
		TS_ASSERT( disulfides.contains( std::make_pair(46,58) ) );
		TS_ASSERT( disulfides.contains( std::make_pair(59,64) ) );
	}

	void test_remove_existing() {
		using namespace protocols::calc_taskop_movers;
		protocols::moves::MoverOP disulf = parse_tag<ForceDisulfidesMover>("<ForceDisulfides name=\"disulf\" remove_existing=\"true\" />");

		core::pose::PoseOP pose = pose_->clone();

		disulf->apply( *pose );

		std::string sequence = pose->annotated_sequence();
		// We're still all set.
		TS_ASSERT_EQUALS(sequence, "L[LEU:NtermProteinFull]TCVTSKSIFGITTENCPDGQNLCFKKWYYIVPRYSDITWGCAATCPKPTNVRETIRCCETDKCNE[GLU:CtermProteinFull]");

		utility::vector1< std::pair<core::Size,core::Size> > disulfides;
		core::conformation::disulfide_bonds( pose->conformation(), disulfides);
		TS_ASSERT_EQUALS( disulfides.size(), 0 );
	}

	void test_rearrange() {
		using namespace protocols::calc_taskop_movers;
		protocols::moves::MoverOP disulf = parse_tag<ForceDisulfidesMover>("<ForceDisulfides name=\"disulf\" disulfides=\"24:64,46A:58,42A:3A\" remove_existing=\"true\" />");

		core::pose::PoseOP pose = pose_->clone();

		disulf->apply( *pose );

		std::string sequence = pose->annotated_sequence();
		// We're still all set.
		TS_ASSERT_EQUALS(sequence, "L[LEU:NtermProteinFull]TC[CYS:disulfide]VTSKSIFGITTENCPDGQNLC[CYS:disulfide]FKKWYYIVPRYSDITWGC[CYS:disulfide]AATC[CYS:disulfide]PKPTNVRETIRC[CYS:disulfide]CETDKC[CYS:disulfide]NE[GLU:CtermProteinFull]");

		utility::vector1< std::pair<core::Size,core::Size> > disulfides;
		core::conformation::disulfide_bonds( pose->conformation(), disulfides);
		TS_ASSERT_EQUALS( disulfides.size(), 3 );
		// We make the somewhat reasonable assumption that the pairs are ordered internally
		// (But we don't assume that any disulfide is in a particular order.
		TS_ASSERT( disulfides.contains( std::make_pair(3,42) ) );
		TS_ASSERT( disulfides.contains( std::make_pair(46,58) ) );
		TS_ASSERT( disulfides.contains( std::make_pair(24,64) ) );
	}

	void test_keep_existing() {
		using namespace protocols::calc_taskop_movers;
		protocols::moves::MoverOP disulf_remove = parse_tag<ForceDisulfidesMover>("<ForceDisulfides name=\"disulf\" disulfides=\"58:46\" remove_existing=\"true\" />");
		protocols::moves::MoverOP disulf = parse_tag<ForceDisulfidesMover>("<ForceDisulfides name=\"disulf\" disulfides=\"3:24,42:17\" />");

		core::pose::PoseOP pose = pose_->clone();

		disulf_remove->apply( *pose );
		utility::vector1< std::pair<core::Size,core::Size> > disulfides;
		core::conformation::disulfide_bonds( pose->conformation(), disulfides);
		TS_ASSERT_EQUALS( disulfides.size(), 1 );
		disulfides.clear();

		disulf->apply( *pose );

		core::conformation::disulfide_bonds( pose->conformation(), disulfides);
		TS_ASSERT_EQUALS( disulfides.size(), 3 );
		// We make the somewhat reasonable assumption that the pairs are ordered internally
		// (But we don't assume that any disulfide is in a particular order.
		TS_ASSERT( disulfides.contains( std::make_pair(46,58) ) );
		TS_ASSERT( disulfides.contains( std::make_pair(3,24) ) );
		TS_ASSERT( disulfides.contains( std::make_pair(17,42) ) );
	}
};
