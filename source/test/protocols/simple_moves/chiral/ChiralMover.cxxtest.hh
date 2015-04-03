// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/chiral/ChiralMover.cxxtest.hh
/// @brief  test suite for protocols::simple_moves::chiral::ChiralMover.cc
/// @author Kevin Drew


// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

// Unit headers
#include <protocols/simple_moves/chiral/ChiralMover.hh>

// Package headers
#include <core/io/pdb/pose_io.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/database/open.hh>

#include <numeric/xyz.functions.hh>

//Auto Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <utility/tools/make_vector1.hh>

static basic::Tracer TR("protocols.simple_moves.chiral.ChiralMover.cxxtest");

// using declarations
using namespace core;
using namespace core::pose;
using namespace protocols;
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace protocols::simple_moves::chiral;
using namespace scoring;

///////////////////////////////////////////////////////////////////////////
/// @name ChiralMoverTest
/// @brief: test switching of residue chirality from L to D
/// @details
///
/// @author Kevin Drew October 17, 2011
///////////////////////////////////////////////////////////////////////////
class ChiralMoverTest : public CxxTest::TestSuite {

public:
	//chemical::ResidueTypeSetCAP residue_set;

	pose::Pose pose;
	ChiralMoverTest() {}

	void setUp() {
		core_init();

		// create score function
		scoring::ScoreFunctionOP score_fxn( ScoreFunctionFactory::create_score_function( scoring::MM_STD_WTS) );

		// read pdb file
		core::import_pose::pose_from_pdb( pose, "protocols/simple_moves/chiral/AAdAdA.pdb" );

		core::init::init_random_generators(1000, "mt19937");
	}

	void tearDown() {

	}

	void test_ChiralMoverLtoL() {

		//core::Real correct_phi ( -131.27 );
		core::Real correct_psi ( 135.00 );

		core::Size lchiral_seq_pos = 1;

		chiral::ChiralMoverOP cm( new chiral::ChiralMover( lchiral_seq_pos, chiral::L_CHIRALITY ) );
		cm->apply( pose );

		//core::Real phi = pose.phi( lchiral_seq_pos );
		core::Real psi = pose.psi( lchiral_seq_pos );
		TR << "psi: " << psi << std::endl;

		// compare to the correct answer
		float TOLERATED_ERROR = 0.1;

		//TS_ASSERT_DELTA( phi, correct_phi, TOLERATED_ERROR );
		TS_ASSERT_DELTA( psi, correct_psi, TOLERATED_ERROR );

		TS_ASSERT( !chiral::is_d_chiral( pose.residue_type( lchiral_seq_pos )));
		TS_ASSERT( chiral::is_l_chiral( pose.residue_type( lchiral_seq_pos )));

		TR << "test_ChiralMoverLtoL completed!! " << std::endl;
	}

	void test_ChiralMoverLtoD() {

		//kdrew: phi and psi should be negated
		core::Real correct_phi ( 139.00 );
		core::Real correct_psi ( -135.00 );

		core::Size lchiral_seq_pos = 2;

		chiral::ChiralMoverOP cm( new chiral::ChiralMover( lchiral_seq_pos, chiral::D_CHIRALITY ) );
		cm->apply( pose );

		core::Real phi = pose.phi( lchiral_seq_pos );
		core::Real psi = pose.psi( lchiral_seq_pos );
		TR << "phi and psi: " << phi << ' ' << psi << std::endl;

		// compare to the correct answer
		float TOLERATED_ERROR = 0.1;

		TS_ASSERT_DELTA( phi, correct_phi, TOLERATED_ERROR );
		TS_ASSERT_DELTA( psi, correct_psi, TOLERATED_ERROR );

		TS_ASSERT( !chiral::is_l_chiral( pose.residue_type( lchiral_seq_pos )));
		TS_ASSERT( chiral::is_d_chiral( pose.residue_type( lchiral_seq_pos )));

		TR << "test_ChiralMoverLtoD completed!! " << std::endl;
	}

	void test_ChiralMoverDtoL() {

		//kdrew: phi and psi should be negated
		core::Real correct_phi ( -139.00 );
		core::Real correct_psi ( 135.00 );

		core::Size dchiral_seq_pos = 3;

		chiral::ChiralMoverOP cm( new chiral::ChiralMover( dchiral_seq_pos, chiral::L_CHIRALITY ) );
		cm->apply( pose );

		core::Real phi = pose.phi( dchiral_seq_pos );
		core::Real psi = pose.psi( dchiral_seq_pos );
		TR << "phi and psi: " << phi << ' ' << psi << std::endl;

		// compare to the correct answer
		float TOLERATED_ERROR = 0.1;

		TS_ASSERT_DELTA( phi, correct_phi, TOLERATED_ERROR );
		TS_ASSERT_DELTA( psi, correct_psi, TOLERATED_ERROR );

		TS_ASSERT( !chiral::is_d_chiral( pose.residue_type( dchiral_seq_pos )));
		TS_ASSERT( chiral::is_l_chiral( pose.residue_type( dchiral_seq_pos )));

		TR << "test_ChiralMoverDtoL completed!! " << std::endl;
	}

	void test_ChiralMoverDtoD() {

		//kdrew: phi and psi should be negated
		core::Real correct_phi ( 139.00 );
		//core::Real correct_psi ( 135.00 );

		core::Size dchiral_seq_pos = 4;

		chiral::ChiralMoverOP cm( new chiral::ChiralMover( dchiral_seq_pos, chiral::D_CHIRALITY ) );
		cm->apply( pose );

		core::Real phi = pose.phi( dchiral_seq_pos );
		//core::Real psi = pose.psi( dchiral_seq_pos );
		TR << "phi : " << phi << std::endl;

		// compare to the correct answer
		float TOLERATED_ERROR = 0.1;

		TS_ASSERT_DELTA( phi, correct_phi, TOLERATED_ERROR );
		//TS_ASSERT_DELTA( psi, correct_psi, TOLERATED_ERROR );

		TS_ASSERT( !chiral::is_l_chiral( pose.residue_type( dchiral_seq_pos )));
		TS_ASSERT( chiral::is_d_chiral( pose.residue_type( dchiral_seq_pos )));

		TR << "test_ChiralMoverDtoD completed!! " << std::endl;
	}

};

