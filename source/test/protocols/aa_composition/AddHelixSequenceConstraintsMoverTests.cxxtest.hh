// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/aa_composition/AddHelixSequenceConstraintsMoverTests.cxxtest.hh
/// @brief  Unit tests for the AddHelixSequenceConstraintsMover.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/aa_composition/AddHelixSequenceConstraintsMover.hh>
#include <protocols/aa_composition/util.hh>

// Protocols Headers
#include <protocols/simple_moves/MutateResidue.hh>

// Core Headers
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/interaction_graph/ResidueArrayAnnealingEvaluator.hh>

// Protocol Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("AddHelixSequenceConstraintsMoverTests");


class AddHelixSequenceConstraintsMoverTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();

		core::pose::PoseOP pose( core::import_pose::pose_from_file( "protocols/aa_composition/2ND2_state1.pdb", false, core::import_pose::PDB_file) );
		for ( core::Size i=1, imax=pose->total_residue(); i <= imax; ++i ) {
			protocols::simple_moves::MutateResidue mut;
			mut.set_target(i);
			mut.set_res_name("GLY");
			mut.apply(*pose);
		}
		testpose_ = pose;
		//pose->dump_pdb("2ND2_state1_glyonly.pdb"); //DELETE ME
	}

	void tearDown(){

	}

	void test_detect_helices() {
		utility::vector1 < std::pair < core::Size, core::Size > > helices;
		protocols::aa_composition::find_helices_over_length( *testpose_, helices, 6 );
		TS_ASSERT_EQUALS( helices.size(), 3 );
		TS_ASSERT_EQUALS( helices[1].first, 4 );
		TS_ASSERT_EQUALS( helices[1].second, 12 );
		TS_ASSERT_EQUALS( helices[2].first, 17 );
		TS_ASSERT_EQUALS( helices[2].second, 29 );
		TS_ASSERT_EQUALS( helices[3].first, 32 );
		TS_ASSERT_EQUALS( helices[3].second, 42 );

		helices.clear();
		//Do it again, with a cutoff that removes one helix.
		protocols::aa_composition::find_helices_over_length( *testpose_, helices, 10 );
		TS_ASSERT_EQUALS( helices.size(), 2 );
		TS_ASSERT_EQUALS( helices[1].first, 17 );
		TS_ASSERT_EQUALS( helices[1].second, 29 );
		TS_ASSERT_EQUALS( helices[2].first, 32 );
		TS_ASSERT_EQUALS( helices[2].second, 42 );

	}

	void test_negative_nterm(){

		core::pose::Pose pose( *testpose_ );

		core::scoring::ScoreFunction scorefxn;
		scorefxn.set_weight( core::scoring::aa_composition, 1 );
		scorefxn.set_weight( core::scoring::ref, 1 );
		utility::vector1< core::Real > ref_wts( core::chemical::num_canonical_aas, 0.0 );
		ref_wts[3] = 1.0; //Give only aspartate a positive reference weight.
		scorefxn.set_method_weights( core::scoring::ref, ref_wts );

		protocols::aa_composition::AddHelixSequenceConstraintsMover add_csts;
		add_csts.set_add_n_terminal_constraints(true, 2, 3, 15.0);
		add_csts.set_add_c_terminal_constraints(false, 2, 3, 15.0);
		add_csts.set_add_overall_constraints(false);
		add_csts.set_add_alanine_constraints(false);
		add_csts.set_add_hydrophobic_constraints(false);
		add_csts.apply(pose);

		(scorefxn)(pose);

		// Set up packer task and packer objects
		core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( pose ));

		utility::vector1< bool > keep_aas( core::chemical::num_canonical_aas, false );
		keep_aas[ core::chemical::aa_ala ] = true;
		keep_aas[ core::chemical::aa_asp ] = true;

		for ( core::Size i = 1; i <= pose.size(); i++ ) {
			task->nonconst_residue_task( i ).restrict_absent_canonical_aas( keep_aas );
		}

		TR << "Pre-design sequence: \t" << pose.sequence() << std::endl;
		TS_ASSERT_EQUALS( pose.sequence(), "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" );

		core::pack::pack_rotamers(pose, scorefxn, task);

		std::string const finalseq( pose.sequence() );
		TR << "Post-design sequence:\t" << pose.sequence() << std::endl;

		core::Size aspcount(0), helix1asp(0), helix2asp(0), helix3asp(0), nontermaspcount(0);
		for ( core::Size i=0, imax( finalseq.length()); i<imax; ++i ) {
			if ( finalseq[i] == 'D' ) {
				++aspcount;
				if ( i >= 3 && i <= 5 ) {
					++helix1asp;
				} else if ( i >= 16 && i <= 18 ) {
					++helix2asp;
				} else if ( i >= 31 && i <= 33 ) {
					++helix3asp;
				} else {
					++nontermaspcount;
				}
			}
		}
		TS_ASSERT_EQUALS( aspcount, 6 );
		TS_ASSERT_EQUALS( helix1asp, 2);
		TS_ASSERT_EQUALS( helix2asp, 2);
		TS_ASSERT_EQUALS( helix3asp, 2);
		TS_ASSERT_EQUALS( nontermaspcount, 0 );
	}

	void test_positive_cterm(){

		core::pose::Pose pose( *testpose_ );

		core::scoring::ScoreFunction scorefxn;
		scorefxn.set_weight( core::scoring::aa_composition, 1 );
		scorefxn.set_weight( core::scoring::ref, 1 );
		utility::vector1< core::Real > ref_wts( core::chemical::num_canonical_aas, 0.0 );
		ref_wts[15] = 1.0; //Give only arginine a positive reference weight.
		scorefxn.set_method_weights( core::scoring::ref, ref_wts );

		protocols::aa_composition::AddHelixSequenceConstraintsMover add_csts;
		add_csts.set_add_n_terminal_constraints(false, 2, 3, 15.0);
		add_csts.set_add_c_terminal_constraints(true, 2, 3, 15.0);
		add_csts.set_add_overall_constraints(false);
		add_csts.set_add_alanine_constraints(false);
		add_csts.set_add_hydrophobic_constraints(false);
		add_csts.apply(pose);

		(scorefxn)(pose);

		// Set up packer task and packer objects
		core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( pose ));

		utility::vector1< bool > keep_aas( core::chemical::num_canonical_aas, false );
		keep_aas[ core::chemical::aa_ala ] = true;
		keep_aas[ core::chemical::aa_arg ] = true;

		for ( core::Size i = 1; i <= pose.size(); i++ ) {
			task->nonconst_residue_task( i ).restrict_absent_canonical_aas( keep_aas );
		}

		TR << "Pre-design sequence: \t" << pose.sequence() << std::endl;
		TS_ASSERT_EQUALS( pose.sequence(), "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" );

		core::pack::pack_rotamers(pose, scorefxn, task);

		std::string const finalseq( pose.sequence() );
		TR << "Post-design sequence:\t" << pose.sequence() << std::endl;

		core::Size argcount(0), helix1arg(0), helix2arg(0), helix3arg(0), nontermargcount(0);
		for ( core::Size i=0, imax( finalseq.length()); i<imax; ++i ) {
			if ( finalseq[i] == 'R' ) {
				++argcount;
				if ( i >= 9 && i <= 11 ) {
					++helix1arg;
				} else if ( i >= 26 && i <= 28 ) {
					++helix2arg;
				} else if ( i >= 39 && i <= 41 ) {
					++helix3arg;
				} else {
					++nontermargcount;
				}
			}
		}
		TS_ASSERT_EQUALS( argcount, 6 );
		TS_ASSERT_EQUALS( helix1arg, 2);
		TS_ASSERT_EQUALS( helix2arg, 2);
		TS_ASSERT_EQUALS( helix3arg, 2);
		TS_ASSERT_EQUALS( nontermargcount, 0 );
	}

	void test_limit_val(){
		core::pose::Pose pose( *testpose_ );

		core::scoring::ScoreFunction scorefxn;
		scorefxn.set_weight( core::scoring::aa_composition, 1 );
		scorefxn.set_weight( core::scoring::ref, 1 );
		utility::vector1< core::Real > ref_wts( core::chemical::num_canonical_aas, 0.0 );
		ref_wts[18] = -1.0; //Give only valine a negative reference weight.
		scorefxn.set_method_weights( core::scoring::ref, ref_wts );

		protocols::aa_composition::AddHelixSequenceConstraintsMover add_csts;
		add_csts.set_add_n_terminal_constraints(false, 2, 3, 15.0);
		add_csts.set_add_c_terminal_constraints(false, 2, 3, 15.0);
		add_csts.set_add_overall_constraints(true, "THR VAL", 2, 15.0);
		add_csts.set_add_alanine_constraints(false);
		add_csts.set_add_hydrophobic_constraints(false);
		add_csts.apply(pose);

		(scorefxn)(pose);

		// Set up packer task and packer objects
		core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( pose ));

		utility::vector1< bool > keep_aas( core::chemical::num_canonical_aas, false );
		keep_aas[ core::chemical::aa_ala ] = true;
		keep_aas[ core::chemical::aa_val ] = true;

		for ( core::Size i = 1; i <= pose.size(); i++ ) {
			task->nonconst_residue_task( i ).restrict_absent_canonical_aas( keep_aas );
		}

		TR << "Pre-design sequence: \t" << pose.sequence() << std::endl;
		TS_ASSERT_EQUALS( pose.sequence(), "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" );

		core::pack::pack_rotamers(pose, scorefxn, task);

		std::string const finalseq( pose.sequence() );
		TR << "Post-design sequence:\t" << pose.sequence() << std::endl;

		core::Size valcount(0), helix1val(0), helix2val(0), helix3val(0), nonhelixvalcount(0);
		for ( core::Size i=0, imax( finalseq.length()); i<imax; ++i ) {
			if ( finalseq[i] == 'V' ) {
				++valcount;
				if ( i >= 3 && i <= 11 ) {
					++helix1val;
				} else if ( i >= 16 && i <= 28 ) {
					++helix2val;
				} else if ( i >= 31 && i <= 41 ) {
					++helix3val;
				} else {
					++nonhelixvalcount;
				}
			}
		}
		TS_ASSERT_EQUALS( valcount, 17 );
		TS_ASSERT_EQUALS( helix1val, 2);
		TS_ASSERT_EQUALS( helix2val, 2);
		TS_ASSERT_EQUALS( helix3val, 2);
		TS_ASSERT_EQUALS( nonhelixvalcount, 11 );

	}

	void test_limit_ala(){
		core::pose::Pose pose( *testpose_ );

		core::scoring::ScoreFunction scorefxn;
		scorefxn.set_weight( core::scoring::aa_composition, 1 );
		scorefxn.set_weight( core::scoring::ref, 1 );
		utility::vector1< core::Real > ref_wts( core::chemical::num_canonical_aas, 0.0 );
		ref_wts[18] = 1.0; //Give only valine a positive reference weight.
		scorefxn.set_method_weights( core::scoring::ref, ref_wts );

		protocols::aa_composition::AddHelixSequenceConstraintsMover add_csts;
		add_csts.set_add_n_terminal_constraints(false, 2, 3, 15.0);
		add_csts.set_add_c_terminal_constraints(false, 2, 3, 15.0);
		add_csts.set_add_overall_constraints(false);
		add_csts.set_add_alanine_constraints(true, 0.5, 100.0, 1.0);
		add_csts.set_add_hydrophobic_constraints(false);
		add_csts.apply(pose);

		(scorefxn)(pose);

		// Set up packer task and packer objects
		core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( pose ));

		utility::vector1< bool > keep_aas( core::chemical::num_canonical_aas, false );
		keep_aas[ core::chemical::aa_ala ] = true;
		keep_aas[ core::chemical::aa_val ] = true;

		for ( core::Size i = 1; i <= pose.size(); i++ ) {
			task->nonconst_residue_task( i ).restrict_absent_canonical_aas( keep_aas );
		}

		TR << "Pre-design sequence: \t" << pose.sequence() << std::endl;
		TS_ASSERT_EQUALS( pose.sequence(), "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" );

		core::pack::pack_rotamers(pose, scorefxn, task);

		std::string const finalseq( pose.sequence() );
		TR << "Post-design sequence:\t" << pose.sequence() << std::endl;

		core::Size helix1val(0), helix2val(0), helix3val(0), nonhelixvalcount(0);
		for ( core::Size i=0, imax( finalseq.length()); i<imax; ++i ) {
			if ( finalseq[i] == 'V' ) {
				if ( i >= 3 && i <= 11 ) {
					++helix1val;
				} else if ( i >= 16 && i <= 28 ) {
					++helix2val;
				} else if ( i >= 31 && i <= 41 ) {
					++helix3val;
				} else {
					++nonhelixvalcount;
				}
			}
		}
		TS_ASSERT_EQUALS( helix1val, 4);
		TS_ASSERT_EQUALS( helix2val, 6);
		TS_ASSERT_EQUALS( helix3val, 5);
		TS_ASSERT_EQUALS( nonhelixvalcount, 0 );

	}

	void test_limit_val_by_hydrophobicity(){
		core::pose::Pose pose( *testpose_ );

		core::scoring::ScoreFunction scorefxn;
		scorefxn.set_weight( core::scoring::aa_composition, 1 );
		scorefxn.set_weight( core::scoring::ref, 1 );
		utility::vector1< core::Real > ref_wts( core::chemical::num_canonical_aas, 0.0 );
		ref_wts[18] = 1.0; //Give only valine a positive reference weight.
		scorefxn.set_method_weights( core::scoring::ref, ref_wts );

		protocols::aa_composition::AddHelixSequenceConstraintsMover add_csts;
		add_csts.set_add_n_terminal_constraints(false, 2, 3, 15.0);
		add_csts.set_add_c_terminal_constraints(false, 2, 3, 15.0);
		add_csts.set_add_overall_constraints(false);
		add_csts.set_add_alanine_constraints(false);
		add_csts.set_add_hydrophobic_constraints(true, 0.5, 1.0);
		add_csts.apply(pose);

		(scorefxn)(pose);

		// Set up packer task and packer objects
		core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( pose ));

		utility::vector1< bool > keep_aas( core::chemical::num_canonical_aas, false );
		keep_aas[ core::chemical::aa_ala ] = true;
		keep_aas[ core::chemical::aa_val ] = true;

		for ( core::Size i = 1; i <= pose.size(); i++ ) {
			task->nonconst_residue_task( i ).restrict_absent_canonical_aas( keep_aas );
		}

		TR << "Pre-design sequence: \t" << pose.sequence() << std::endl;
		TS_ASSERT_EQUALS( pose.sequence(), "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" );

		core::pack::pack_rotamers(pose, scorefxn, task);

		std::string const finalseq( pose.sequence() );
		TR << "Post-design sequence:\t" << pose.sequence() << std::endl;

		core::Size helix1val(0), helix2val(0), helix3val(0), nonhelixvalcount(0);
		for ( core::Size i=0, imax( finalseq.length()); i<imax; ++i ) {
			if ( finalseq[i] == 'V' ) {
				if ( i >= 3 && i <= 11 ) {
					++helix1val;
				} else if ( i >= 16 && i <= 28 ) {
					++helix2val;
				} else if ( i >= 31 && i <= 41 ) {
					++helix3val;
				} else {
					++nonhelixvalcount;
				}
			}
		}
		TS_ASSERT_EQUALS( helix1val, 5 );
		TS_ASSERT_EQUALS( helix2val, 7 );
		TS_ASSERT_EQUALS( helix3val, 6 );
		TS_ASSERT_EQUALS( nonhelixvalcount, 0 );

	}


private:

	core::pose::PoseCOP testpose_;

};



