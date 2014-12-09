// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/optimization/Minimizer.cxxtest.hh
/// @brief  test suite for Minimizer
/// @author Phil Bradley
/// @author Sergey Lyskov

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>

// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

// AUTO-REMOVED #include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <core/id/AtomID_Mask.hh>
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.optimization.Minimizer.cxxtest");

using namespace core;

class MinimizerTests : public CxxTest::TestSuite
{
	chemical::ResidueTypeSetCAP residue_set;

public:
	MinimizerTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();

		residue_set = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );
	}

	// Shared finalization goes here.
	void tearDown() {
	}


///////////////////////////////////////////////////////////////////////////////
// ------------------------------------------ //
/// @brief simple test minimization
void test_simple_min()
{
	using namespace optimization;
	using pose::Pose;
	using id::AtomID;
	using id::DOF_ID;
	using id::PHI;
	using id::THETA;
	using id::D;

	pose::Pose start_pose(create_test_in_pdb_pose());
	//core::import_pose::pose_from_pdb( start_pose, "core/optimization/test_in.pdb" );

	kinematics::MoveMapOP mm( new kinematics::MoveMap );

	// setup moving dofs
	for ( int i=30; i<= 35; ++i ) {
		mm->set_bb ( i, true );
		mm->set_chi( i, true );
	}

	// setup the options
	scoring::ScoreFunctionOP scorefxn( new scoring::ScoreFunction );

	AtomTreeMinimizer minimizer;
	MinimizerOptionsOP min_options( new MinimizerOptions( "linmin", 10.0, true, true, false ) );

// non-core level code has no place in a "core" test!
// furthermore, core does not (and should not) link the protocols lib!
// move this test to protocols if you want to include/test code from protocols!!
//	protocols::simple_moves::MinMover min_mover( mm, scorefxn, "linmin", 10.0, true
//		/*use_nblist*/, true /*deriv_check*/, false /*no verbose-deriv-check, default*/ );

	{ // test out a couple different scoring functions

		{ // just fa_rama
			scorefxn->set_weight( scoring::rama, 1.0 );

			Pose pose;
			pose = start_pose;
			TR << "MINTEST: rama" << std::endl;
			minimizer.run( pose, *mm, *scorefxn, *min_options );
		}

		/*{ // just fa_elec // apl -- temporarily disable fa_elec
			scorefxn->reset();
			scorefxn->set_weight( scoring::fa_elec, 0.5 );

			Pose pose;
			pose = start_pose;
			TR << "MINTEST: fa_elec" << std::endl;
			minimizer.run( pose, *mm, *scorefxn, *min_options );
		}*/

		{ // just fa_dun
			scorefxn->reset();
			scorefxn->set_weight( scoring::fa_dun, 1.0 );

			Pose pose;
			pose = start_pose;
			TR << "MINTEST: fa_dun" << std::endl;
			minimizer.run( pose, *mm, *scorefxn, *min_options );
		}

		{ // just fa_atr
			scorefxn->reset();
			scorefxn->set_weight( scoring::fa_atr, 0.80 );

			Pose pose;
			pose = start_pose;
			TR << "MINTEST: atr" << std::endl;
			minimizer.run( pose, *mm, *scorefxn, *min_options );

		}

		{ // just fa_atr, rep, sol
			scorefxn->reset();
			scorefxn->set_weight( scoring::fa_atr, 0.80 );
			scorefxn->set_weight( scoring::fa_rep, 0.44 );
			scorefxn->set_weight( scoring::fa_sol, 0.65 );

			Pose pose;
			pose = start_pose;
			TR << "MINTEST: atr-rep-sol" << std::endl;
			minimizer.run( pose, *mm, *scorefxn, *min_options );
		}

		{ // fa_atr, rep, sol and fa_intra atr, rep, & sol
			scorefxn->reset();
			scorefxn->set_weight( scoring::fa_atr, 0.80 );
			scorefxn->set_weight( scoring::fa_rep, 0.44 );
			scorefxn->set_weight( scoring::fa_sol, 0.65 );

			scorefxn->set_weight( scoring::fa_intra_atr, 0.80 );
			scorefxn->set_weight( scoring::fa_intra_rep, 0.44 );
			scorefxn->set_weight( scoring::fa_intra_sol, 0.65 );

			Pose pose;
			pose = start_pose;
			TR << "MINTEST: atr-rep-sol and intra atr-rep-sol" << std::endl;
			TR << "start score: " << (*scorefxn)( pose ) << std::endl;
			minimizer.run( pose, *mm, *scorefxn, *min_options );
			pose.dump_pdb( "min_intrares.pdb" );
			TR << "end score: " << (*scorefxn)( pose ) << std::endl;
		}

		{  // p_aa_pp
			scorefxn->reset();
			scorefxn->set_weight( scoring::p_aa_pp, 0.29 );

			Pose pose;
			pose = start_pose;
			TR << "MINTEST: p_aa_pp" << std::endl;
			TR << "start score: " << (*scorefxn)( pose ) << std::endl;
			minimizer.run( pose, *mm, *scorefxn, *min_options );
			pose.dump_pdb( "min_intrares.pdb" );
			TR << "end score: " << (*scorefxn)( pose ) << std::endl;
		}

		{ // fa_atr, rep, sol and p_aa_pp
			scorefxn->reset();
			scorefxn->set_weight( scoring::fa_atr, 0.80 );
			scorefxn->set_weight( scoring::fa_rep, 0.44 );
			scorefxn->set_weight( scoring::fa_sol, 0.65 );
			scorefxn->set_weight( scoring::p_aa_pp, 0.29 );

			Pose pose;
			pose = start_pose;
			TR << "MINTEST: atr-rep-sol and p_aa_pp" << std::endl;
			TR << "start score: " << (*scorefxn)( pose ) << std::endl;
			minimizer.run( pose, *mm, *scorefxn, *min_options );
			pose.dump_pdb( "min_intrares.pdb" );
			TR << "end score: " << (*scorefxn)( pose ) << std::endl;
		}


		{ // just fa_atr, rep, sol, rigid-body minimization
			scorefxn->reset();
			scorefxn->set_weight( scoring::fa_atr, 0.80 );
			scorefxn->set_weight( scoring::fa_rep, 0.44 );
			scorefxn->set_weight( scoring::fa_sol, 0.65 );

			Pose pose;
			pose = start_pose;

			{ // setup a foldtree
				kinematics::FoldTree f( pose.total_residue() );
				f.new_jump( 8, 26, 18 );
				f.reorder( 8 );
				pose.fold_tree( f );
			}

			kinematics::MoveMapOP mm2( new kinematics::MoveMap );
			mm2->set_jump( 1, true );

			TR << "MINTEST: atr-rep-sol jumpmin" << std::endl;
			minimizer.run( pose, *mm2, *scorefxn, *min_options );

		}

		{ // just fa_pair
			scorefxn->reset();
			scorefxn->set_weight( scoring::fa_pair, 1.0 );

			Pose pose;
			pose = start_pose;
			TR << "MINTEST: fa_pair" << std::endl;
			minimizer.run( pose, *mm, *scorefxn, *min_options );
		}

		{ // just backbone hbonds
			scorefxn->reset();
			scorefxn->set_weight( scoring::hbond_lr_bb, 1.0 );
			scorefxn->set_weight( scoring::hbond_sr_bb, 1.0 );

			Pose pose;
			pose = start_pose;
			TR << "MINTEST: bb hbonds" << std::endl;
			minimizer.run( pose, *mm, *scorefxn, *min_options );
		}

		{ // all hbonds
			scorefxn->reset();
			scorefxn->set_weight( scoring::hbond_lr_bb, 1.0 );
			scorefxn->set_weight( scoring::hbond_sr_bb, 1.0 );
			scorefxn->set_weight( scoring::hbond_bb_sc, 1.0 );
			scorefxn->set_weight( scoring::hbond_sc, 1.0 );

			Pose pose;
			pose = start_pose;
			TR << "MINTEST: all hbonds" << std::endl;
			minimizer.run( pose, *mm, *scorefxn, *min_options );
		}

		return;

	} // scope


	Pose pose;
	pose = start_pose;

	// set the moving dofs
	kinematics::MoveMapOP mm1( new kinematics::MoveMap );
	kinematics::MoveMapOP mm2( new kinematics::MoveMap );
	kinematics::MoveMapOP mm3( new kinematics::MoveMap );
	kinematics::MoveMapOP mm4( new kinematics::MoveMap );
	kinematics::MoveMapOP mm5( new kinematics::MoveMap );
	// single backbone
	mm1->set_bb( 4, true );

	// all bb and chi
	mm2->set_bb( true );
	mm2->set_chi( true );

	// single dof
	mm3->set( DOF_ID( AtomID(1,4), PHI ), true );

	// everything!
	mm4->set( PHI, true );
	mm4->set( THETA, true );
	mm4->set( D, true );


	// everything! + neighborlist auto-update
	mm5->set( PHI, true );
	mm5->set( THETA, true );
	mm5->set( D, true );

	// setup scorefxn
	scorefxn->reset();
	scorefxn->set_weight( scoring::fa_atr, 0.80 );

	pose.dump_pdb( "before.pdb" );
	MinimizerOptionsOP min_options2( new MinimizerOptions( "dfpmin", 0.001, true, true, false ) );
	minimizer.run( pose, *mm1, *scorefxn, *min_options2 );

	pose.dump_pdb( "after1.pdb" );

	exit(0);

	minimizer.run( pose, *mm2, *scorefxn, *min_options2 );
	pose.dump_pdb( "after2.pdb" );

	minimizer.run( pose, *mm3, *scorefxn, *min_options2 );
	pose.dump_pdb( "after3.pdb" );

	minimizer.run( pose, *mm4, *scorefxn, *min_options2 );
	pose.dump_pdb( "after4.pdb" );

	min_options2->nblist_auto_update( true );
	minimizer.run( pose, *mm5, *scorefxn, *min_options2 );
	pose.dump_pdb( "after5.pdb" );

};


};
