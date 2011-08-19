// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/util/cout.cxxtest.hh
/// @brief  testing presense of working 'cout' statements in the code
/// @author Sergey Lyskov

// Test headers
#include <cxxtest/TestSuite.h>

#include "platform/types.hh"

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <core/scoring/ScoreFunction.hh>
#include <protocols/moves/ScoreMover.hh>
#include <protocols/moves/BackboneMover.hh>

// Package Headers

#include <core/io/pdb/pose_io.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/NamedStubID.hh>
#include <core/io/pdb/file_data.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/signals/DestructionEvent.hh>
#include <core/pose/signals/EnergyEvent.hh>
#include <core/scoring/types.hh>
#include <core/scoring/constraints/Constraints.fwd.hh>
#include <utility/stream_util.hh>
#include <protocols/moves/MoverStatistics.hh>
#include <utility/keys/Key2Tuple.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/format.hh>


static basic::Tracer TR("protocol.util.cout.cxxtest");

using namespace core;

class cout_io : public CxxTest::TestSuite
{
	chemical::ResidueTypeSetCAP residue_set;

public:
	cout_io() {}


	// Shared initialization goes here.
	void setUp() {
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// ---------------------------------------------------------------------------------------------
	/// @brief testing presence of uncommented  cout statements in the code.
	///        If this test failed that mean that some of the testes functions have working 'cout'
	///        inside them - replace it with Tracer output to fix the test.
	void test_cout_io() {
		basic::Tracer::super_mute(true);  // Mute all channels no mater what...

		std::stringbuf sb;
		std::streambuf * cout_buff = std::cout.rdbuf();
		// intersepting cout streambuf
		std::cout.rdbuf(&sb);


		// testig Pose IO
		core_init_with_additional_options( "-no_optH" );
		residue_set = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );
		core::pose::Pose pose;
		core::import_pose::pose_from_pdb(pose, "protocols/util/test_in.pdb");
		core::io::pdb::dump_pdb(pose, "._tmp_pdb_.pdb");
		core::pose::Pose P2(pose);

		scoring::ScoreFunction scorefxn;
		scorefxn.set_weight(scoring::fa_atr, 1.0 );
		scorefxn( pose );

		{ protocols::moves::ScoreMover sm;  sm.test_move(pose); }
		{ protocols::moves::SmallMover sm;  sm.test_move(pose); }
		{ protocols::moves::ShearMover sm;  sm.test_move(pose); }


		// setting back old cout streambuf
		std::cout.rdbuf(cout_buff);
		//std::cout << "sb=" << sb.str() << std::endl;
		basic::Tracer::super_mute(false);  // unmute all channels

		if( sb.str() != "" ) {
			std::cout << "Pose code contain cout statements: " << sb.str() << std::endl;
			TS_FAIL("");
		}
	}
};

