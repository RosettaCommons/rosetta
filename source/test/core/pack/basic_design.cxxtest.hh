// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/basic_design.cxxtest.hh
/// @brief  test suite for resfile reader
/// @author Steven Lewis

// Test headers
#include <cxxtest/TestSuite.h>

#include <core/types.hh>

#include <test/core/init_util.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/PackerTask_.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/palette/DefaultPackerPalette.hh>
#include <core/pack/palette/CustomBaseTypePackerPalette.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pack/task/RotamerSampleOptions.hh>

#include <core/pack/task/ResfileReader.hh>
#include <core/pack/pack_rotamers.hh>

#include <string>
#include <sstream> //stringstreams can convert integers into strings type-safely for comparisons en masse

//Auto Headers
#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR("core.pack.basic_design.cxxtest");

//I'm lazy using's
using namespace core;
using namespace pack;
using namespace task;
using namespace palette;
using namespace pose;
using namespace chemical;
using namespace utility::pointer;
using std::string;
using std::stringstream;

// --------------- Test Class --------------- //

class BasicDesignTests : public CxxTest::TestSuite {

public:

	Pose pose;

	// --------------- Fixtures --------------- //

	void setUp() {
		core_init();
		core::pose::make_pose_from_sequence( pose, "AAA", "fa_standard" );
	}

	void tearDown() {
		pose.clear(); //nuke that sucker in case it got altered
		//nothing necessary for OP packer task - setUp's regeneration of a fresh copy will cause the old OP's data to delete itself


	}

	// --------------- Test Cases --------------- //

	void test_design_rotamers(){
		core::scoring::ScoreFunction sfxn;

		DefaultPackerPaletteOP pp = make_shared< DefaultPackerPalette >();
		PackerTask_OP the_task = make_shared< PackerTask_ >( pose, pp ); // ALSO has an implicit DPP
		the_task->show(TR);
		TS_ASSERT( the_task->residue_task( 1 ).allowed_residue_types().size() == 21 );
	}

	void test_design_custom(){
		core::scoring::ScoreFunction sfxn;

		CustomBaseTypePackerPaletteOP pp = make_shared< CustomBaseTypePackerPalette >();
		pp->add_type( "NVL" );
		PackerTask_OP the_task = make_shared< PackerTask_ >( pose, pp ); // ALSO has an implicit DPP
		the_task->show(TR);
		TS_ASSERT( the_task->residue_task( 1 ).allowed_residue_types().size() == 22 );
	}

	void test_design_rotamers_TF(){
		core::scoring::ScoreFunction sfxn;

		//the_task->initialize_from_command_line();
		DefaultPackerPaletteOP pp = make_shared< DefaultPackerPalette >();
		TaskFactory tf;
		tf.set_packer_palette( pp );
		PackerTaskOP the_task = tf.create_task_and_apply_taskoperations( pose );
		the_task->show(TR);
		TS_ASSERT( the_task->residue_task( 1 ).allowed_residue_types().size() == 21 );
	}

	void test_design_custom_TF(){
		core::scoring::ScoreFunction sfxn;

		CustomBaseTypePackerPaletteOP pp = make_shared< CustomBaseTypePackerPalette >();
		pp->add_type( "NVL" );

		TaskFactory tf;
		tf.set_packer_palette( pp );
		PackerTaskOP the_task = tf.create_task_and_apply_taskoperations( pose );
		the_task->show(TR);
		TS_ASSERT( the_task->residue_task( 1 ).allowed_residue_types().size() == 22 );
	}

};//end class
