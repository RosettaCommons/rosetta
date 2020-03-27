// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_moves/FoldTreeFromMotif.cxxtest.hh
/// @brief  test suite for FoldTreeFromLoops
/// @author Brian Coventry (bcov@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>

#include <platform/types.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <protocols/rosetta_scripts/XmlObjects.hh>
#include <protocols/simple_moves/FoldTreeFromMotif.hh>


#include <core/types.hh>


#include <test/UTracer.hh>
#include <test/core/init_util.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>

using namespace core;


static basic::Tracer TR("protocols.simple_moves.FoldTreeFromMotif.cxxtest");

class FoldTreeFromMotifTests : public CxxTest::TestSuite
{

public:
	void setUp()
	{
		core_init();
	}

	void test_fold_tree_from_loops() {

		core::pose::PoseOP pose_p = utility::pointer::make_shared<core::pose::Pose>();
		import_pose::pose_from_file( *pose_p, "protocols/motif_grafting/2nm1_1zx3a_hybrid.pdb" );

		core::pose::Pose pose = *(pose_p->clone());


		protocols::simple_moves::FoldTreeFromMotif ftfm;
		core::select::residue_selector::ResidueSelectorOP sel = utility::pointer::make_shared<core::select::residue_selector::ResidueIndexSelector>("5,6,7");

		ftfm.residue_selector(sel);

		ftfm.apply(pose);

		TS_ASSERT( pose.fold_tree().root() == 6 );

		protocols::rosetta_scripts::XmlObjectsCOP objs = protocols::rosetta_scripts::XmlObjects::create_from_string(
			"<RESIDUE_SELECTORS><Index name=\"sel\" resnums=\"5,6,7\"/></RESIDUE_SELECTORS>\n"
			"<MOVERS><FoldTreeFromMotif name=\"ftfm2\" residue_selector=\"sel\" /></MOVERS>"
		);

		protocols::moves::MoverOP ftfm2 = objs->get_mover("ftfm2");

		pose = *(pose_p->clone());

		ftfm2->apply(pose);

		TS_ASSERT( pose.fold_tree().root() == 6 );

	}

};
