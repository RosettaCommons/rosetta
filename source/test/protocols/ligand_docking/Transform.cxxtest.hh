// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/conformation/UltraLightResidue.hh>
#include <protocols/qsar/scoring_grid/GridManager.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <protocols/ligand_docking/Transform.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>
//#include <test/core/init_util.hh>


static basic::Tracer TR("protocols.ligand_docking.Transform.cxxtest");


class TransformTests : public CxxTest::TestSuite {

private:

	core::pose::Pose pose_;
	protocols::ligand_docking::Transform mover_;
	core::Size begin_;

public:

	void setUp() {
		core_init_with_additional_options("-extra_res_fa protocols/ligand_docking/ZNx.params protocols/ligand_docking/7cpa.params");

		core::import_pose::pose_from_file( pose_, "protocols/ligand_docking/7cpa_7cpa_native.pdb" , core::import_pose::PDB_file);
		begin_ = pose_.conformation().chain_begin(core::pose::get_chain_id_from_chain('X', pose_));

	}

	void tearDown() {}

	void test_score_constraints(){
		using namespace protocols::ligand_docking;
		core::conformation::UltraLightResidue test_ligand(pose_.residue(begin_).get_self_ptr());

		core::pose::Pose copy_pose = pose_;
		core::scoring::func::HarmonicFuncOP func( new core::scoring::func::HarmonicFunc( 2.5 , 0.5 ) );

		//  core::Size begin=pose_.conformation().chain_begin(core::pose::get_chain_id_from_chain('X', pose_));

		//Actual distance between identified atoms is 2.262014
		core::scoring::constraints::AtomPairConstraintOP cst1( new core::scoring::constraints::AtomPairConstraint(
			core::id::AtomID( pose_.residue( 248 ).atom_index( "OH" ), 248 ),
			core::id::AtomID( pose_.residue( begin_ ).atom_index( "H32"   ), begin_ ),
			func ) );

		copy_pose.add_constraint( cst1 );

		core::scoring::ScoreFunctionOP sfxn(new core::scoring::ScoreFunction);
		sfxn->set_weight( core::scoring::atom_pair_constraint, 1.0 );

		TS_ASSERT_DELTA(mover_.score_constraints(copy_pose, test_ligand, sfxn), 0.226549, 0.001);
	}

	void test_randomize_ligand() {
		using namespace protocols::ligand_docking;

		core::conformation::UltraLightResidue test_ligand(pose_.residue(begin_).get_self_ptr());
		core::conformation::UltraLightResidue start_ligand(test_ligand);

		numeric::random::rg().set_seed( "mt19937", time(0) );

		core::Vector start_center(start_ligand.center());

		//Test initial perturb of residue

		for ( core::Size i=0; i <= 100; i++ ) {
			mover_.randomize_ligand(test_ligand, 5, 360);
			TS_ASSERT_LESS_THAN(test_ligand.center().distance(start_center), 5);
			test_ligand = start_ligand;
		}

	}

	void test_conformer_change(){

		//Test conformer change of residue

		core::conformation::UltraLightResidue test_ligand(pose_.residue(begin_).get_self_ptr());
		core::conformation::UltraLightResidue start_ligand(test_ligand);

		core::Vector start_center(start_ligand.center());

		core::Size begin=pose_.conformation().chain_begin(core::pose::get_chain_id_from_chain('X', pose_));

		mover_.setup_conformers(pose_, begin);

		for ( core::Size i=0; i <= 100; i++ ) {
			mover_.change_conformer(test_ligand);
			TS_ASSERT_LESS_THAN(test_ligand.center().distance(start_center), 5);
			test_ligand = start_ligand;
		}


	}

};

