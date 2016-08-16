// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
// libRosetta headers

#include <devel/cartesian_frags/CartesianFragment.hh>
#include <devel/cartesian_frags/dna_util.hh>
#include <devel/cartesian_frags/DNA_FragLib.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Conformation.hh>

#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularPowerFunc.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>
#include <basic/prof.hh> // profiling
#include <basic/Tracer.hh>

#include <numeric/constants.hh>

#include <utility/vector1.hh>
#include <numeric/xyz.functions.hh>


// // // C++ headers

// //silly using/typedef
namespace devel {
namespace cartesian_frags {

using namespace core;
using utility::vector1;

/// tracer object:
static THREAD_LOCAL basic::Tracer tt( "devel.cartesian_frags.dna_util", basic::t_trace );

///////////////////////////////////////////////////////////////////////////////

/// @details  Optimize the bond angles and torsions of a sugar->sugar backbone suite to match a desired transform
///   uses gradient based minimization with a tether to starting DOF values.

void
optimize_suite(
	Direction const & dir,
	kinematics::RT const & target,
	pose::Pose & pose, // the mini pose
	CartesianFragment & frag
)
{
	using namespace pose;
	using namespace scoring;
	using namespace conformation;
	using namespace chemical;
	using namespace scoring::constraints;
	using namespace kinematics;
	using namespace id;
	using namespace optimization;

	// tethering parameters:
	Real const torsion_constraint_weight( 1.0 );
	Real const torsion_constraint_sdev_degrees( 10.0 );
	int  const torsion_constraint_power( 8 );
	Real const angle_constraint_weight( 5.0 );
	Real const angle_constraint_sdev_degrees( 2.0 );
	int  const angle_constraint_power( 2 );

	Real const coordinate_constraint_weight( 1.0 );
	Real const        dof_constraint_weight( 0.1 );

	assert( pose.total_residue() == 4 );

	RT const & rt( frag.rt(1) );

	Stub const start_stub( rt );
	Stub const target_stub( target );

	//// setup frag atom ids //////////////////////////////////
	vector1< SafeAtomID > frag_atoms;
	frag_atoms.push_back( SafeAtomID( "O3*", 0 ) );
	frag_atoms.push_back( SafeAtomID( "P"  , 1 ) );
	frag_atoms.push_back( SafeAtomID( "O5*", 1 ) );
	frag_atoms.push_back( SafeAtomID( "C5*", 1 ) );
	if ( dir == Backward ) {
		std::reverse( frag_atoms.begin(), frag_atoms.end() );
		for ( int i=1; i<= 4; ++i ) --frag_atoms[i].seqpos;
	}

	frag_atoms.insert( frag_atoms.begin(), SafeAtomID( "instub_atom2", 0 ) );
	frag_atoms.insert( frag_atoms.begin(), SafeAtomID( "instub_atom3", 0 ) );

	frag_atoms.push_back( SafeAtomID( "outstub1_atom2", 0 ) );
	frag_atoms.push_back( SafeAtomID( "outstub1_atom3", 0 ) );
	///// setup my atom ids ///////////////////////////////////////
	vector1< AtomID > my_atoms;
	for ( Size i=1; i<=3; ++i ) my_atoms.push_back( AtomID( i, 1   ) );
	for ( Size i=1; i<=2; ++i ) my_atoms.push_back( AtomID( 1, i+1 ) );
	for ( Size i=1; i<=3; ++i ) my_atoms.push_back( AtomID( i, 4   ) );
	//// sanity checks
	assert( frag.xyz( frag_atoms[3] ).length()    < 1e-3 ); // 1st atom should be at origin
	// last atom should be at rt translation (outgoing stub origin)
	assert( frag.xyz( frag_atoms[6] ).distance(rt.get_translation()) < 1e-3 );
	//// copy frag coords to pose:
	for ( Size i=1; i<= my_atoms.size(); ++i ) pose.set_xyz( my_atoms[i], frag.xyz( frag_atoms[ i ] ) );
	pose::Pose start_pose;
	start_pose = pose;

	// setup the options
	MinimizerOptions options( "lbfgs_armijo_nonmonotone", 0.001, true /*use_nblist*/, false /*deriv_check*/,
		false /*no verbose-deriv-check, is default*/ );
	//// setup the moving degrees of freedom
	kinematics::MoveMap mm;

	vector1< DOF_ID > angles, torsions;

	torsions.push_back( DOF_ID( AtomID( 1, 2 ), PHI   ) );
	torsions.push_back( DOF_ID( AtomID( 1, 3 ), PHI   ) );
	torsions.push_back( DOF_ID( AtomID( 1, 4 ), PHI   ) );
	torsions.push_back( DOF_ID( AtomID( 2, 4 ), PHI   ) );
	torsions.push_back( DOF_ID( AtomID( 3, 4 ), PHI   ) );

	angles.push_back( DOF_ID( AtomID( 1, 2 ), THETA ) );
	angles.push_back( DOF_ID( AtomID( 1, 3 ), THETA ) );
	angles.push_back( DOF_ID( AtomID( 1, 4 ), THETA ) );
	angles.push_back( DOF_ID( AtomID( 2, 4 ), THETA ) );

	for ( Size i=1; i<= torsions.size(); ++i ) mm.set( torsions[i], true );
	for ( Size i=1; i<=   angles.size(); ++i ) mm.set(   angles[i], true );
	//// setup constraints for outgoing stub + to tether angles to starting values
	ConstraintSetOP cst_set( new ConstraintSet() );
	core::scoring::func::FuncOP coord_cst_func( new core::scoring::func::HarmonicFunc( 0.0, 0.1 ) );
	for ( Size i=1; i<= 3; ++i ) {
		cst_set->add_constraint( ConstraintCOP( ConstraintOP( new CoordinateConstraint( AtomID(i,4), AtomID(1,1), target_stub.build_fake_xyz(i),
			coord_cst_func ) ) ) );
	}

	for ( Size i=1; i<= torsions.size(); ++i ) {
		Real const target( pose.dof( torsions[i] ) );
		Real const sdev( numeric::conversions::radians( torsion_constraint_sdev_degrees ) );
		Real const weight( torsion_constraint_weight );
		cst_set->add_dof_constraint( torsions[i], core::scoring::func::FuncOP( new core::scoring::func::CircularPowerFunc( target, sdev, torsion_constraint_power, weight ) ));
	}

	for ( Size i=1; i<= angles.size(); ++i ) {
		Real const target( pose.dof( angles[i] ) );
		Real const sdev( numeric::conversions::radians( angle_constraint_sdev_degrees ) );
		Real const weight( angle_constraint_weight );
		cst_set->add_dof_constraint( angles[i], core::scoring::func::FuncOP( new core::scoring::func::CircularPowerFunc( target, sdev, angle_constraint_power, weight ) ));
	}

	//// setup scorefxn
	ScoreFunction scorefxn;
	scorefxn.set_weight( coordinate_constraint, coordinate_constraint_weight );
	scorefxn.set_weight( dof_constraint, dof_constraint_weight );

	pose.constraint_set( cst_set );

	//// run minimizer!!
	AtomTreeMinimizer().run( pose, mm, scorefxn, options );

	//// now update the coords in frag, the outgoing frag

	for ( Size i=1; i<= torsions.size(); ++i ) {
		frag.set_torsion_angle( frag_atoms[i], frag_atoms[i+1], frag_atoms[i+2], frag_atoms[i+3],
			pose.dof( torsions[i] ) );
	}

	for ( Size i=1; i<= angles.size(); ++i ) {
		frag.set_bond_angle( frag_atoms[i+1], frag_atoms[i+2], frag_atoms[i+3],
			numeric::constants::d::pi - pose.dof( angles[i] ) );
	}

	for ( Size i=1; i<= my_atoms.size(); ++i ) {
		assert( pose.xyz( my_atoms[i] ).distance( frag.xyz( frag_atoms[i] ) ) < 1e-2 );
	}

}
///////////////////////////////////////////////////////////////////////////////
/// @details  Scan a fragment library for fragment combinations that patch up the DNA backbone.
///   Specifically: given target transforms from the DNA base forward and backward to the next sugar,
///   find sugar+fwd-suite+bwd-suite combos that optimally close.
///   Useful if you've moved the base and want to fix any chainbreaks. Calls optimize_suite to
///   perfectly close the best fragment combination.

Real
find_sugar_and_suite_frags(
	DNA_FragLib const & lib,
	kinematics::RT const & fwd_target,
	kinematics::RT const & bwd_target,
	bool const lower_terminus,
	bool const upper_terminus,
	CartesianFragment & sugar,
	CartesianFragment & fwd_suite,
	CartesianFragment & bwd_suite
)
{
	using kinematics::RT;

	Real best_summed_dev( 999.9 );

	PROF_START( basic::FIND_SUGAR_AND_SUITE_FRAGS_I );
	for ( Size ii=1; ii<= lib.sugars.size(); ++ii ) {
		RT const & s1( lib.sugars[ii].stub_transform(2) ); //  forward to epsilon
		RT const & s2( lib.sugars[ii].stub_transform(1) ); // backward to gamma

		// get the target suite transforms
		RT const fwd_suite_target( s1, fwd_target );
		RT const bwd_suite_target( s2, bwd_target );

		vector1< Size > best_index( 2, 0 );
		vector1< Real > best_dev( 2, 999.9 );

		assert( lib.forward_suites.size() == lib.backward_suites.size() );
		for ( Size jj=1; jj<= lib.forward_suites.size(); ++jj ) {
			for ( int r=1; r<=2; ++r ) {
				if ( ( upper_terminus && r == 1 ) || ( lower_terminus && r == 2 ) ) continue;
				CartesianFragment const & frag( r == 1 ? lib.forward_suites[jj] : lib.backward_suites[jj] );
				RT const & target( r == 1 ? fwd_suite_target : bwd_suite_target );

				Real const dev( frag.rt( 1 ).distance_squared( target ) );
				if ( dev < best_dev[r] ) {
					best_dev[r] = dev;
					best_index[r] = jj;
				}
			}
		}

		if ( upper_terminus ) { best_dev[1] = 0.0; best_index[1] = 1; }
		if ( lower_terminus ) { best_dev[2] = 0.0; best_index[2] = 1; }

		if ( best_dev[1] + best_dev[2] < best_summed_dev ) {
			best_summed_dev = best_dev[1] + best_dev[2];
			sugar = lib.sugars[ii];
			fwd_suite =  lib.forward_suites[ best_index[1] ];
			bwd_suite = lib.backward_suites[ best_index[2] ];
		}
	}
	PROF_STOP( basic::FIND_SUGAR_AND_SUITE_FRAGS_I );

	PROF_START( basic::FIND_SUGAR_AND_SUITE_FRAGS_II );
	Real final_summed_dev( 0.0 );
	{ // try optimizing the suites

		RT const & s1( sugar.stub_transform(2) ); //  forward to epsilon
		RT const & s2( sugar.stub_transform(1) ); // backward to gamma

		// get the target suite transforms
		RT const fwd_suite_target( s1, fwd_target );
		RT const bwd_suite_target( s2, bwd_target );

		if ( !upper_terminus ) {
			optimize_suite(  Forward, fwd_suite_target, lib.suite_pose(), fwd_suite );
			final_summed_dev += fwd_suite_target.distance_squared( fwd_suite.rt(1) );
		}
		if ( !lower_terminus ) {
			optimize_suite( Backward, bwd_suite_target, lib.suite_pose(), bwd_suite );
			final_summed_dev += bwd_suite_target.distance_squared( bwd_suite.rt(1) );
		}
	}
	PROF_STOP( basic::FIND_SUGAR_AND_SUITE_FRAGS_II );
	tt << "find_sugar_and_suite_frags:: final_summed_dev: " << final_summed_dev << std::endl;
	return final_summed_dev;
}
///////////////////////////////////////////////////////////////////////////////
/// @details  Patches up the backbone before and/or after DNA position i by varying the sugar and both suites
//
//
Real
patch_up_backbone(
	Size const i,
	DNA_FragLib const & lib,
	conformation::Conformation & conf
)
{
	using namespace kinematics;
	using namespace id;
	Stub chi1, epsilon, epsilon_p, gamma, gamma_n;

	bool const lower_terminus( conf.residue(i).is_lower_terminus() );
	bool const upper_terminus( conf.residue(i).is_upper_terminus() );

	chi1    = torsion_stub( TorsionID( i  , CHI, 1 ), Backward, conf );
	epsilon = torsion_stub( TorsionID( i  ,  BB, 5 ),  Forward, conf );
	gamma   = torsion_stub( TorsionID( i  ,  BB, 3 ), Backward, conf );
	if ( !lower_terminus ) epsilon_p = torsion_stub( TorsionID( i-1,  BB, 5 ),  Forward, conf );
	if ( !upper_terminus ) gamma_n   = torsion_stub( TorsionID( i+1,  BB, 3 ), Backward, conf );

	RT const target_forward ( chi1, gamma_n   ); // forward to next sugar
	RT const target_backward( chi1, epsilon_p ); // backward to previous sugar

	CartesianFragment sugar, fwd_suite, bwd_suite;
	Real const final_summed_dev
		( find_sugar_and_suite_frags( lib, target_forward, target_backward, lower_terminus, upper_terminus,
		sugar, fwd_suite, bwd_suite ) );

	// now make the fragment insertions
	sugar.insert( conf, i );
	if ( !upper_terminus ) fwd_suite.insert( conf, i );
	if ( !lower_terminus ) bwd_suite.insert( conf, i );

	return final_summed_dev;

}

///////////////////////////////////////////////////////////////////////////////
/// @details  Patches up the DNA backbone link between i and i+1
//
void
patch_up_backbone_link(
	Size const i,
	DNA_FragLib const & lib,
	scoring::ScoreFunction const &, // scorefxn,
	pose::Pose & pose_inout
)
{

	//pose::Pose pose( pose_inout );

	if ( numeric::random::uniform() < 0.5 ) {
		// try i then i+1
		patch_up_backbone( i  , lib, pose_inout.conformation() );
		patch_up_backbone( i+1, lib, pose_inout.conformation() );
	} else {
		// try i+1 then i
		patch_up_backbone( i+1, lib, pose_inout.conformation() );
		patch_up_backbone( i  , lib, pose_inout.conformation() );
	}
	//  Real const fwd_score( scorefxn( pose ) );
	//  Real const bwd_score( scorefxn( pose_inout ) );
	//  tt << "patch_up_backbone_link: " << i << ' ' << fwd_score << ' ' << bwd_score << std::endl;
	//  if ( fwd_score < bwd_score ) {
	//   pose_inout = pose;
	//  }
}
} // ns cartesian_frags
} // ns devel
