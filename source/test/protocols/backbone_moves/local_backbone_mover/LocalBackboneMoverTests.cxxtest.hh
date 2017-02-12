// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/backbone_moves/local_backbone_mover/LocalBackboneMoverTests.cxxtest.hh
/// @brief  Test the BackboneMover and its subordinate classes
/// @author xingjiepan (xingjiepan@gmail.com)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/backbone_moves/local_backbone_mover/FreePeptide.hh>
#include <protocols/backbone_moves/local_backbone_mover/GapCloser.hh>
#include <protocols/backbone_moves/local_backbone_mover/free_peptide_movers/TranslationFreePeptideMover.hh>
#include <protocols/backbone_moves/local_backbone_mover/free_peptide_movers/LongAxisRotationFreePeptideMover.hh>
#include <protocols/backbone_moves/local_backbone_mover/free_peptide_movers/ShearFreePeptideMover.hh>
#include <protocols/backbone_moves/local_backbone_mover/free_peptide_movers/CircularPermuteFreePeptideMover.hh>

// Numeric Headers
#include <numeric/xyzVector.io.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>

// Core Headers
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Protocol Headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("LocalBackboneMoverTests");
using namespace protocols::backbone_moves::local_backbone_mover;


class LocalBackboneMoverTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();
		pose_op_ = core::import_pose::pose_from_file("protocols/backbone_moves/local_backbone_mover/1arb.pdb");
	}

	void tearDown(){
	}

	void test_free_peptide(){
		using core::id::AtomID;
		
		//pose_op_->dump_file("1arb.pdb");
		FreePeptide fp(*pose_op_, 185, 189);

		TS_ASSERT_DELTA( fp.n_ca_bond(187), 1.5115, 1E-3);
		TS_ASSERT_DELTA( fp.ca_c_bond(187), 1.5023, 1E-3);
		TS_ASSERT_DELTA( fp.c_n_bond(187), 1.3156, 1E-3);
		TS_ASSERT_DELTA( fp.n_ca_c_angle(187), 1.9197, 1E-3);
		TS_ASSERT_DELTA( fp.ca_c_n_angle(187), 2.0423, 1E-3);
		TS_ASSERT_DELTA( fp.c_n_ca_angle(187), 2.0260, 1E-3);
		TS_ASSERT_DELTA( fp.phi(187), -72.0884, 1E-3); 
		TS_ASSERT_DELTA( fp.psi(187), 158.0593, 1E-3); 
		TS_ASSERT_DELTA( fp.omega(187), -179.1394, 1E-3); 

		for(core::Size i=187; i < 189; ++i){
			for(core::Size j=1; j < 4; ++j){
				AtomID aid(j, i);

				numeric::xyzVector <core::Real> pos1 = pose_op_->xyz(aid);
				numeric::xyzVector <core::Real> pos2;
				switch(j){
					case 1:
						pos2 = fp.n_xyz(i);
						break;
					case 2:
						pos2 = fp.ca_xyz(i);
						break;
					case 3:
						pos2 = fp.c_xyz(i);
						break;
				}

				for(core::Size k=1; k < 4; ++k){
					TS_ASSERT_DELTA( pos1(k), pos2(k), 1E-3 );
				}
			}
		}

		fp.n_ca_bond(187, 3);
		//fp.phi(187, 180);

		fp.apply_to_pose(*pose_op_);
		//pose_op_->dump_file("1arb_after_apply_free_peptide.pdb");
	}

	void test_free_peptide_move(){
	
		FreePeptide fp(*pose_op_, 185, 189);
		
		// Test translate

		xyzVector <Real> t(0, 2, 0);
		fp.translate(t);
		TS_ASSERT_DELTA(fp.ca_xyz(187)(2), -1.485, 0.01);

		// Test rotate

		fp.rotate(numeric::x_rotation_matrix(3.14159));
		TS_ASSERT_DELTA(fp.ca_xyz(187)(2), 0.979, 0.01);
	
		// Test align
		
		fp.align();	
		TS_ASSERT_DELTA(fp.ca_xyz(187)(2), -3.485, 0.01);
	}

	void test_gap_closer(){
		FreePeptide fp(*pose_op_, 185, 189);
		GapCloser gc;

		gc.solve_gaps(fp);

		TS_ASSERT(gc.gap_solved());

		Real expected_torsions1[2][6] = {{298.566, 333.011, 233.702, 19.9289, 273.951, 184.753},
			{321.719, 295.247, 267.547, 314.165, 359.09, 173.648}};

		for(Size i=1; i <=gc.pivot_torsions1_.size(); ++i){
			for(Size j=1; j <= 6; ++j){
				TS_ASSERT_DELTA(gc.pivot_torsions1_[i][j], expected_torsions1[i - 1][j - 1], 1e-2);
			}
		}
	
		Real expected_torsions2[2][6] = {{232.145, 210.428, 224.168, 159.234, 221.898, 170.051},
			{236.955, 236.523, 213.899, 237.808, 145.54, 154.649}};
		
		for(Size i=1; i <=gc.pivot_torsions2_.size(); ++i){
			for(Size j=1; j <= 6; ++j){
				TS_ASSERT_DELTA(gc.pivot_torsions2_[i][j], expected_torsions2[i - 1][j - 1], 1e-2);
			}
		}

		gc.apply_closure(*pose_op_, fp);
		//pose_op_->dump_file("1arb_after_close_gaps.pdb");
	}

	void test_free_peptide_movers(){
		FreePeptide fp(*pose_op_, 185, 189);
		GapCloser gc;

		free_peptide_movers::FreePeptideMoverOP fpm(
			//new free_peptide_movers::TranslationFreePeptideMover(0.5));
			//new free_peptide_movers::LongAxisRotationFreePeptideMover(1, true));
			//new free_peptide_movers::ShearFreePeptideMover(186, 60));
			new free_peptide_movers::CircularPermuteFreePeptideMover(1, true));

		fpm->apply(fp);

		gc.solve_gaps(fp);
		TS_ASSERT(gc.gap_solved());
		gc.apply_closure(*pose_op_, fp);
		
		//pose_op_->dump_file("1arb_after_bb_move.pdb");
	}
private:
	core::pose::PoseOP pose_op_;

};



