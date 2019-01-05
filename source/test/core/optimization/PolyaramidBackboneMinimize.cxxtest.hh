// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/Minimizer.cxxtest.hh
/// @brief  test suite for Minimizer
/// @author Phil Bradley
/// @author Sergey Lyskov

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>


#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Mask.hh>
#include <utility/vector1.hh>


using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.optimization.PolyaramidBackboneMinimizeTests.cxxtest");

using namespace core;
using namespace core::pose;
using namespace core::optimization;
using namespace core::scoring;
using namespace core::scoring::func;
using namespace core::scoring::constraints;
using namespace core::id;

class PolyaramidBackboneMinimizeTests : public CxxTest::TestSuite
{
	chemical::ResidueTypeSetCOP residue_set;
	Pose start_pose;
	ScoreFunctionOP sfxn = nullptr;

public:
	PolyaramidBackboneMinimizeTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();

		residue_set = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );

		sfxn = ScoreFunctionOP( ScoreFunctionFactory::create_score_function( "empty.wts" ) );
		sfxn->set_weight( dihedral_constraint, 1.0 );
	}

	Pose init_pose( std::string const & resn ) {
		start_pose = Pose();
		start_pose.append_residue_by_jump(core::conformation::Residue( residue_set->name_map( resn ), true ), 1 );
		start_pose.append_residue_by_bond(core::conformation::Residue( residue_set->name_map( resn ), true ), true );
		start_pose.append_residue_by_bond(core::conformation::Residue( residue_set->name_map( resn ), true ), true );

		for ( Size ii = 1; ii <= 3; ++ii ) {
			start_pose.set_phi( ii, 20 );
			start_pose.set_psi( ii, 20 );
			start_pose.set_omega( ii, 180 );
		}
		return start_pose;
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	/// @brief Can we minimize a tightly constrained polyaramid phi?
	/// @author Andy Watkins (amw579@stanford.edu)
	void test_min_phi_succeeds_ortho_ala() {

		Pose test_pose = init_pose( "ORTHO_POLYARAMID_ALA" );

		// constrain phi
		test_pose.add_constraint( DihedralConstraintOP( new DihedralConstraint(
			AtomID( test_pose.residue_type( 1 ).atom_index( "C" ), 1),
			AtomID( test_pose.residue_type( 2 ).atom_index( "N" ), 2),
			AtomID( test_pose.residue_type( 2 ).atom_index( "CB1" ), 2),
			AtomID( test_pose.residue_type( 2 ).atom_index( "CA" ), 2),
			CircularHarmonicFuncOP( new CircularHarmonicFunc( 0, 0.01 ) ) ) ) );

		TR << "score start: " << ( *sfxn )( test_pose ) << std::endl;

		kinematics::MoveMapOP mm( new kinematics::MoveMap );

		// Who cares about anything but the relevant dihedral really
		mm->set_bb( 2, 1, true );

		AtomTreeMinimizer minimizer;
		MinimizerOptionsOP min_options( new MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.01, true, false, false ) );

		minimizer.run( test_pose, *mm, *sfxn, *min_options );

		TR << "score end: " << ( *sfxn )( test_pose ) << std::endl;

		test_pose.dump_pdb( "done_min.pdb" );
		TS_ASSERT_DELTA( test_pose.phi( 2 ), 0, 1e-4 );
	}

	/// @brief Can we minimize a tightly constrained polyaramid psi?
	/// @author Andy Watkins (amw579@stanford.edu)
	void test_min_psi_succeeds_ortho_ala() {

		Pose test_pose = init_pose( "ORTHO_POLYARAMID_ALA" );

		// constrain psi
		test_pose.add_constraint( DihedralConstraintOP( new DihedralConstraint(
			AtomID( test_pose.residue_type( 2 ).atom_index( "CB1" ), 2),
			AtomID( test_pose.residue_type( 2 ).atom_index( "CA" ), 2),
			AtomID( test_pose.residue_type( 2 ).atom_index( "C" ), 2),
			AtomID( test_pose.residue_type( 3 ).atom_index( "N" ), 3),
			CircularHarmonicFuncOP( new CircularHarmonicFunc( 0, 0.01 ) ) ) ) );

		TR << "score start: " << ( *sfxn )( test_pose ) << std::endl;

		kinematics::MoveMapOP mm( new kinematics::MoveMap );

		// Who cares about anything but the relevant dihedral really
		mm->set_bb( 2, 3, true );

		AtomTreeMinimizer minimizer;
		MinimizerOptionsOP min_options( new MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.01, true, false, false ) );

		minimizer.run( test_pose, *mm, *sfxn, *min_options );

		TR << "score end: " << ( *sfxn )( test_pose ) << std::endl;

		test_pose.dump_pdb( "done_min.pdb" );
		TS_ASSERT_DELTA( test_pose.psi( 2 ), 0, 1e-4 );
	}

	/// @brief Can we minimize a tightly constrained polyaramid phi?
	/// @author Andy Watkins (amw579@stanford.edu)
	void test_min_phi_succeeds_meta_ala() {

		Pose test_pose = init_pose( "META_POLYARAMID_ALA" );

		// constrain phi
		test_pose.add_constraint( DihedralConstraintOP( new DihedralConstraint(
			AtomID( test_pose.residue_type( 1 ).atom_index( "C" ), 1),
			AtomID( test_pose.residue_type( 2 ).atom_index( "N" ), 2),
			AtomID( test_pose.residue_type( 2 ).atom_index( "CG1" ), 2),
			AtomID( test_pose.residue_type( 2 ).atom_index( "CB1" ), 2),
			CircularHarmonicFuncOP( new CircularHarmonicFunc( 0, 0.01 ) ) ) ) );

		TR << "score start: " << ( *sfxn )( test_pose ) << std::endl;

		kinematics::MoveMapOP mm( new kinematics::MoveMap );

		// Who cares about anything but the relevant dihedral really
		mm->set_bb( 2, 1, true );

		AtomTreeMinimizer minimizer;
		MinimizerOptionsOP min_options( new MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.01, true, false, false ) );

		minimizer.run( test_pose, *mm, *sfxn, *min_options );

		TR << "score end: " << ( *sfxn )( test_pose ) << std::endl;

		test_pose.dump_pdb( "done_min.pdb" );
		TS_ASSERT_DELTA( test_pose.phi( 2 ), 0, 1e-4 );
	}

	/// @brief Can we minimize a tightly constrained polyaramid psi?
	/// @author Andy Watkins (amw579@stanford.edu)
	void test_min_psi_succeeds_meta_ala() {

		Pose test_pose = init_pose( "META_POLYARAMID_ALA" );

		// constrain psi
		test_pose.add_constraint( DihedralConstraintOP( new DihedralConstraint(
			AtomID( test_pose.residue_type( 2 ).atom_index( "CB1" ), 2),
			AtomID( test_pose.residue_type( 2 ).atom_index( "CA" ), 2),
			AtomID( test_pose.residue_type( 2 ).atom_index( "C" ), 2),
			AtomID( test_pose.residue_type( 3 ).atom_index( "N" ), 3),
			CircularHarmonicFuncOP( new CircularHarmonicFunc( 0, 0.01 ) ) ) ) );

		TR << "score start: " << ( *sfxn )( test_pose ) << std::endl;

		kinematics::MoveMapOP mm( new kinematics::MoveMap );

		// Who cares about anything but the relevant dihedral really
		mm->set_bb( 2, 4, true );

		AtomTreeMinimizer minimizer;
		MinimizerOptionsOP min_options( new MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.01, true, false, false ) );

		minimizer.run( test_pose, *mm, *sfxn, *min_options );

		TR << "score end: " << ( *sfxn )( test_pose ) << std::endl;

		test_pose.dump_pdb( "done_min.pdb" );
		TS_ASSERT_DELTA( test_pose.psi( 2 ), 0, 1e-4 );
	}

	/// @brief Can we minimize a tightly constrained polyaramid phi?
	/// @author Andy Watkins (amw579@stanford.edu)
	void test_min_phi_succeeds_para_ala() {

		Pose test_pose = init_pose( "PARA_POLYARAMID_ALA" );

		// constrain phi
		test_pose.add_constraint( DihedralConstraintOP( new DihedralConstraint(
			AtomID( test_pose.residue_type( 1 ).atom_index( "C" ), 1),
			AtomID( test_pose.residue_type( 2 ).atom_index( "N" ), 2),
			AtomID( test_pose.residue_type( 2 ).atom_index( "CD" ), 2),
			AtomID( test_pose.residue_type( 2 ).atom_index( "CG1" ), 2),
			CircularHarmonicFuncOP( new CircularHarmonicFunc( 0, 0.01 ) ) ) ) );

		TR << "score start: " << ( *sfxn )( test_pose ) << std::endl;

		kinematics::MoveMapOP mm( new kinematics::MoveMap );

		// Who cares about anything but the relevant dihedral really
		mm->set_bb( 2, 1, true );

		AtomTreeMinimizer minimizer;
		MinimizerOptionsOP min_options( new MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.01, true, false, false ) );

		minimizer.run( test_pose, *mm, *sfxn, *min_options );

		TR << "score end: " << ( *sfxn )( test_pose ) << std::endl;

		test_pose.dump_pdb( "done_min.pdb" );
		TS_ASSERT_DELTA( test_pose.phi( 2 ), 0, 1e-4 );
	}

	/// @brief Can we minimize a tightly constrained polyaramid psi?
	/// @author Andy Watkins (amw579@stanford.edu)
	void test_min_psi_succeeds_para_ala() {

		Pose test_pose = init_pose( "PARA_POLYARAMID_ALA" );

		// constrain psi
		test_pose.add_constraint( DihedralConstraintOP( new DihedralConstraint(
			AtomID( test_pose.residue_type( 2 ).atom_index( "CB1" ), 2),
			AtomID( test_pose.residue_type( 2 ).atom_index( "CA" ), 2),
			AtomID( test_pose.residue_type( 2 ).atom_index( "C" ), 2),
			AtomID( test_pose.residue_type( 3 ).atom_index( "N" ), 3),
			CircularHarmonicFuncOP( new CircularHarmonicFunc( 0, 0.01 ) ) ) ) );

		TR << "score start: " << ( *sfxn )( test_pose ) << std::endl;

		kinematics::MoveMapOP mm( new kinematics::MoveMap );

		// Who cares about anything but the relevant dihedral really
		mm->set_bb( 2, 5, true );

		AtomTreeMinimizer minimizer;
		MinimizerOptionsOP min_options( new MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.01, true, false, false ) );

		minimizer.run( test_pose, *mm, *sfxn, *min_options );

		TR << "score end: " << ( *sfxn )( test_pose ) << std::endl;

		test_pose.dump_pdb( "done_min.pdb" );
		TS_ASSERT_DELTA( test_pose.psi( 2 ), 0, 1e-4 );
	}

};
