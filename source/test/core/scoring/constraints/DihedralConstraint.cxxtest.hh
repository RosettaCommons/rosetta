// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constaints/DihedralConstraint.cxxtest.hh
/// @brief  test suite for dihedral constraints
/// @author Ian Davis

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>
#include <test/util/deriv_funcs.hh>

#include <core/conformation/Residue.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>


#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/FourPointsFunc.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/func/XYZ_Func.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

#include <numeric/conversions.hh>

//Auto Headers
#include <core/id/AtomID_Mask.hh>
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.constraints.DihedralConstraint.cxxtest");

using namespace core;

class DihedralConstraintTests : public CxxTest::TestSuite
{

public:
	DihedralConstraintTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	void test_dihe_score() {
		using namespace core;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		core::scoring::func::FourPointsFunc fourpts;
		fourpts.xyz( 1, Vector( 0, 0, 0 ) );
		fourpts.xyz( 2, Vector( 0, 1.0, 0 ));
		fourpts.xyz( 3, Vector( 0.707, 0.707, 0 ));
		fourpts.xyz( 4, Vector( 0.707, 0.707, 1.0 )); // 90 degrees

		core::scoring::func::CircularHarmonicFuncOP func = new core::scoring::func::CircularHarmonicFunc( numeric::conversions::radians( 80 ), 10 );

		AtomID at1( 1, 1), at2( 2, 1 ), at3( 3, 1 ), at4( 4, 1 );

		DihedralConstraint dcst( at1, at2, at3, at4, func );
		EnergyMap weights, emap;
		weights[ dihedral_constraint ] = 1.0;
		dcst.score( fourpts, weights, emap );
		//Size before_precision = std::cout.precision();
		//std::cout.precision( 16 );
		//std::cout << "Dihedral constraint func: " << emap[ dihedral_constraint ] << std::endl;
		TS_ASSERT_DELTA( emap[ dihedral_constraint ],  0.02467401100272339, 1e-12 );

		//std::cout.precision( before_precision );
	}

	void test_dihedral_derivatives() {
		using namespace core;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		core::pose::PoseOP ubqstump = create_twores_1ubq_poseop();
		TS_ASSERT( ubqstump->total_residue() == 2 );
		AtomID at1( 1, 1), at2( 2, 1 ), at3( 3, 1 ), at4( 4, 1 );
		core::scoring::func::CircularHarmonicFuncOP func = new core::scoring::func::CircularHarmonicFunc( numeric::conversions::radians( 80 ), 10 );

		ScoreFunction sfxn;
		sfxn.set_weight( dihedral_constraint, 1.0 );
		kinematics::MoveMap movemap;
    movemap.set_bb( true );
    movemap.set_chi( true );

		for ( Size ii = 1; ii <= ubqstump->total_residue(); ++ii ) {
			core::chemical::ResidueType const & rsd_type( ubqstump->residue_type( ii ));
			// for each dihedral angle in the residue type
			for ( Size dihe = 1; dihe <= rsd_type.ndihe(); ++dihe ) {
				ubqstump->remove_constraints();

				at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = ii;
				at1.atomno() = ( rsd_type.dihedral( dihe ) ).key1();
				at2.atomno() = ( rsd_type.dihedral( dihe ) ).key2();
				at3.atomno() = ( rsd_type.dihedral( dihe ) ).key3();
				at4.atomno() = ( rsd_type.dihedral( dihe ) ).key4();

				DihedralConstraintOP dcst = new DihedralConstraint( at1, at2, at3, at4, func );
				ubqstump->add_constraint( dcst );

				AtomDerivValidator adv( *ubqstump, sfxn, movemap );
				/// This call runs a numeric deriv check on all the free dofs in the system and makes sure
				/// that the analytic norm matches the numeric norm to 1e-3.
				adv.simple_deriv_check( true, 1e-6 );
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////
	// ------------------------------------------ //
	/// @brief simple test minimization
	void test_dihe_constraints()
	{
		using namespace core::scoring::constraints;
		using core::chemical::ResidueType;
		using core::chemical::AtomIndices;
		using core::conformation::Residue;
		using core::id::AtomID;

		pose::Pose start_pose( create_test_in_pdb_pose() );
		//core::import_pose::pose_from_pdb( start_pose, "core/scoring/constraints/test_in.pdb" );
		pose::Pose pose = start_pose; // a copy

		scoring::ScoreFunctionOP scorefxn = new scoring::ScoreFunction;
		scorefxn->reset();
		scorefxn->set_weight( scoring::fa_atr, 0.80 );
		scorefxn->set_weight( scoring::fa_rep, 0.44 );
		scorefxn->set_weight( scoring::fa_sol, 0.65 );
		scorefxn->set_weight( scoring::dihedral_constraint, 1.0 );

		// input stddev is in degrees, but dihedral constraints deal in radians
		core::Real const stddev_radians = numeric::conversions::radians( 5.0 );
		ConstraintSetOP constraints = new ConstraintSet();

		core::Size const which_res = 116; // Leu 164
		core::Size const which_chi = 2;
		Residue const & rsd = pose.residue(which_res);
		ResidueType const & rsd_type = pose.residue_type(which_res);
		TR << "Constraint on rsd " << rsd.name3() << " " << which_res << ", chi " << which_chi << std::endl;

		core::Real const start_chi_degrees = rsd.chi(which_chi);
		core::Real const start_chi_radians = numeric::conversions::radians( start_chi_degrees );
		core::scoring::func::FuncOP restr_func = new core::scoring::func::CircularHarmonicFunc( start_chi_radians, stddev_radians );
		AtomIndices chi_idx = rsd_type.chi_atoms(which_chi); // 1-based
		ConstraintOP constraint = new DihedralConstraint(
			AtomID(chi_idx[1], which_res),
			AtomID(chi_idx[2], which_res),
			AtomID(chi_idx[3], which_res),
			AtomID(chi_idx[4], which_res),
			restr_func
		);
		constraints->add_constraint( constraint );
		pose.constraint_set( constraints );
		TR << "Constraint: " << constraint->atom(1) << " " << constraint->atom(2) << " " << constraint->atom(3) << " " << constraint->atom(4) << std::endl;
		TR << "Starting position " << rsd.chi(which_chi) << " == " << start_chi_degrees << " degrees" << std::endl;

		(*scorefxn)( pose );
		TR << "Starting constraint energy " << pose.energies().total_energies()[ scoring::dihedral_constraint ] << " (expected: 0)" << std::endl;
		TS_ASSERT_DELTA( 0.0, pose.energies().total_energies()[ scoring::dihedral_constraint ], 1e-6 );

		//rsd.chi( rsd.chi()+20.0 );
		pose.set_chi( which_chi, which_res, rsd.chi(which_chi)+20.0 );
		TR << "New position " << rsd.chi(which_chi) << " degrees" << std::endl;
		(*scorefxn)( pose );
		TR << "New constraint energy " << pose.energies().total_energies()[ scoring::dihedral_constraint ] << " (expected: >0)" << std::endl;
		TS_ASSERT( 10.0 < pose.energies().total_energies()[ scoring::dihedral_constraint ] );

		kinematics::MoveMapOP mm = new kinematics::MoveMap;
		mm->set_chi( which_res, true );

		TR << "Minimizing..." << std::endl;
	// BAD!  We are in core, not protocols!  protocols lib is not even linked here!
	//protocols::simple_moves::MinMover min_mover( mm, scorefxn, "dfpmin_armijo_nonmonotone_atol", 0.001, true /*use_nblist*/ );

		core::optimization::AtomTreeMinimizer minimizer;
		core::optimization::MinimizerOptionsOP min_options =
			new core::optimization::MinimizerOptions( "dfpmin_armijo_nonmonotone_atol", 0.001, true, true, false );
		minimizer.run( pose, *mm, *scorefxn, *min_options );

	//	min_mover.apply( pose );
		(*scorefxn)( pose );
		TR << "New position " << rsd.chi(which_chi) << " degrees" << std::endl;
		TS_ASSERT_DELTA( start_chi_degrees, rsd.chi(which_chi), 0.5 );
		TR << "New constraint energy " << pose.energies().total_energies()[ scoring::dihedral_constraint ] << " (expected: ~0)" << std::endl;
		TS_ASSERT_DELTA( 0.0, pose.energies().total_energies()[ scoring::dihedral_constraint ], 1e-1 );

		//[old behavior] When a pose is copied, we make a (very) deep copy of its constraints.
		// When a pose is copied, we make a shallow copy of its constraint set
		// Check that all constraints are copied and compare as equal.
		pose::Pose pose_copy( pose );
		TS_ASSERT( pose.constraint_set() != pose_copy.constraint_set() );
		utility::vector1< ConstraintCOP > constraints_orig = pose.constraint_set()->get_all_constraints();
		utility::vector1< ConstraintCOP > constraints_copy = pose_copy.constraint_set()->get_all_constraints();
		TS_ASSERT( constraints_orig.size() == 1 );
		TS_ASSERT( constraints_orig.size() == constraints_copy.size() );
		for(core::Size i = 1; i <= constraints_orig.size(); ++i) {
			TR << "Original " << i << ": " << constraints_orig[i]->to_string() << std::endl;
			TR << "Copy " << i << ": " << constraints_copy[i]->to_string() << std::endl;
			//[old behavior] Objects are equivalent (same data) but not identical (different places in memory)
			//TS_ASSERT( constraints_orig[i] != constraints_copy[i] );
			// Objects are identical (same place in memory)
			TS_ASSERT( constraints_orig[i] == constraints_copy[i] );
			TS_ASSERT( constraints_orig[i]() == constraints_copy[i]() );
			TS_ASSERT( constraints_orig[i]->to_string() == "DihedralConstraint 116 2 116 6 116 7 116 8 CircularHarmonicFunc 2.44839 0.0872665" );
			TS_ASSERT( constraints_orig[i]->to_string() == constraints_copy[i]->to_string() );
		}

		core::Size const size_before = pose.constraint_set()->get_all_constraints().size();
		TS_ASSERT( size_before == 1 );
		TS_ASSERT( pose.remove_constraint(constraint) == true );
		core::Size const size_after = pose.constraint_set()->get_all_constraints().size();
		TS_ASSERT( size_after == 0 );
		// Can't remove a constraint twice:
		TS_ASSERT( pose.remove_constraint(constraint) == false );
	}


};
