// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/loop_graph/evaluator/SixDTransRotPotentialEvaluator.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/scoring/loop_graph/evaluator/SixDTransRotPotentialEvaluator.hh>
#include <core/scoring/loop_graph/util.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/constants.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "core.scoring.loop_graph.evaluator.SixDTransRotPotentialEvaluator" );

using namespace core;
using numeric::constants::r::pi;

namespace core {
namespace scoring {
namespace loop_graph {
namespace evaluator {

using Matrix = numeric::xyzMatrix<Real>;

//Constructor
SixDTransRotPotentialEvaluator::SixDTransRotPotentialEvaluator( Size const & takeoff_pos,
	Size const & landing_pos,
	pose::Pose const & pose,
	core::Real const & loop_fixed_cost,
	SixDTransRotPotential const & potential ):
	LoopClosePotentialEvaluator( log( 1.0e27 / 6.022e23) + loop_fixed_cost ), // note factor to convert from 1M reference state of 6D potentials to units of A^-3.
	takeoff_pos_( takeoff_pos ),
	landing_pos_( landing_pos ),
	potential_( potential )
{
	runtime_assert( ( landing_pos_ - takeoff_pos_ ) > 0 ); // will need to update if/when we handle cyclization.
	set_loop_closure_energy( evaluate( pose ) );
	figure_out_if_loop_involves_current_pose( pose );
}

//Destructor
SixDTransRotPotentialEvaluator::~SixDTransRotPotentialEvaluator() = default;

///////////////////////////////////////////////////////////////////////////////////////////////
/// @details currently limited to RNA loops but would be easy to generalize to proteins, etc.
Real
SixDTransRotPotentialEvaluator::evaluate( pose::Pose const & pose )
{
	using namespace core::id;
	using namespace core::kinematics;

	AtomID takeoff_a_atom_id, takeoff_b_atom_id, takeoff_c_atom_id;
	AtomID landing_a_atom_id, landing_b_atom_id, landing_c_atom_id;
	Vector takeoff_xyz, takeoff_a_xyz, takeoff_b_xyz, takeoff_c_xyz;
	Vector landing_xyz, landing_a_xyz, landing_b_xyz, landing_c_xyz;

	// see, e.g., RNA_FragmentMonteCarlo output_jump_o3p_to_o5p mode for saving jump information.
	get_loop_atom( takeoff_pos_, pose, " O3'", takeoff_atom_id_,  takeoff_xyz );
	get_loop_atom( takeoff_pos_, pose, " O3'", takeoff_a_atom_id, takeoff_a_xyz );
	get_loop_atom( takeoff_pos_, pose, " C3'", takeoff_b_atom_id, takeoff_b_xyz );
	get_loop_atom( takeoff_pos_, pose, " C4'", takeoff_c_atom_id, takeoff_c_xyz );

	get_loop_atom( landing_pos_, pose, " O5'", landing_atom_id_,  landing_xyz );
	get_loop_atom( landing_pos_, pose, " C5'", landing_a_atom_id, landing_a_xyz );
	get_loop_atom( landing_pos_, pose, " O5'", landing_b_atom_id, landing_b_xyz );
	get_loop_atom( landing_pos_, pose, " C4'", landing_c_atom_id, landing_c_xyz );

	// takeoff
	Stub stub1( takeoff_xyz /* center */,
		takeoff_a_xyz /* a */,
		takeoff_b_xyz /* b  [b->a defines x] */,
		takeoff_c_xyz /* c  [c->b defines y] */ );
	stub1.M = Matrix::cols( stub1.M.col_y(), stub1.M.col_z(), stub1.M.col_x() ); // Prefer to have C3'->O3' (takeoff vector) along z

	// landing
	Stub stub2( landing_xyz /* center */,
		landing_a_xyz /* a */,
		landing_b_xyz /* b  [b->a defines x] */,
		landing_c_xyz /* c  [c->b defines y] */ );
	stub2.M = Matrix::cols( stub2.M.col_y(), stub2.M.col_z(), stub2.M.col_x() ); // Prefer to have O5'->C5' (landing vector) along z

	stub1_ = stub1;
	j_ = Jump( stub1, stub2 );

	Real const score = loop_fixed_cost() + potential_.evaluate( j_ );
	return score;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
void
SixDTransRotPotentialEvaluator::figure_out_if_loop_involves_current_pose( pose::Pose const & pose )
{
	using namespace core::pose::full_model_info;
	set_involves_current_pose( false );
	// for derivs, need to flag whether moving an AtomID in the 'current pose' will affect the potential.
	if ( full_model_info_defined( pose ) ) {
		if ( const_full_model_info( pose ).res_list().has_value( takeoff_pos_ ) ) {
			runtime_assert( const_full_model_info( pose ).res_list().has_value( landing_pos_ ) );
			set_involves_current_pose( true );
		}
	} else { // no full_model_info.
		set_involves_current_pose( true );
	}
	if ( involves_current_pose() ) {
		set_current_pose_takeoff_atom( takeoff_atom_id_ );
		set_current_pose_landing_atom( landing_atom_id_ );
		set_current_pose_takeoff_xyz( pose.xyz( takeoff_atom_id_ ) );
		set_current_pose_landing_xyz( pose.xyz( landing_atom_id_ ) );
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// yes, I actually worked this out using quaternions & stuff, and checked in MATLAB.
//
// See:
//
//  https://github.com/rhiju/loop_close/blob/master/notes/Das_6DPotential_Derivs.pdf
//  https://github.com/rhiju/loop_close/blob/master/notes/Das_6DPotential_F1_F2.pdf
//  https://github.com/rhiju/loop_close/tree/master/check_rotations
//
// and check accompanying repo.
//
// also see comments in code below.
//
// -- rhiju, 2017
void
SixDTransRotPotentialEvaluator::get_f1_f2( Vector & f1, Vector & f2, bool const takeoff ) const
{
	std::pair< Vector, Vector > const deriv = potential_.get_derivative( j_ );
	Vector const & dEdr = deriv.first;

	/////////////////////////////////////////////////////////////////////
	// First, contributions from translation derivative.
	//
	// This is the 'standard' form of the f1 and f2, where
	//
	//   f2 = dE/dr   and  f1 = f2 x r
	//
	// as used elsewhere in the code, from Abe and Go.
	// [these vectors are equivalent to force & torque in classical
	//  mechanics, up to a minus sign.]
	//
	// E is defined with respect to the transform going from
	//  the takeoff coordinate system to the landing coordinate system
	//
	// So, straight out of the potential,  dE/dr is in the local frame of the
	//  takeoff coordinate system. We need to rotate it back to
	//  the global coordinate frame:
	/////////////////////////////////////////////////////////////////////
	Vector const dEdr_rot = stub1_.rotation() * dEdr;
	f2 = -dEdr_rot;
	// The cross product is with the landing atom position -- pretty
	//  easy to check this if you imagine a torsion angle change in
	//  the chain that intervenes between takeoff and landing, a la
	//  Abe and Go:
	f1 = cross( f2, current_pose_landing_xyz() );

	/////////////////////////////////////////////////////////////////////
	// second, contributions from rotational function.
	/////////////////////////////////////////////////////////////////////
	Vector const dEdv = deriv.second * ( 180.0 / pi ); // convert to radians
	Vector const rotvector = numeric::rotation_axis_angle( j_.get_rotation() ); // in radians
	Real const V( rotvector.length() );
	Vector f1_rot( 0.0 );
	if ( V > 0.0 ) {
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// This transforms the rotation derivative dE/dv (where v is the rotation vector)
		//  to the appropriate F1 ('torque') to use in the global frame.
		//
		// Derived based on chain rule and quaternion rotations (see link above
		//  to scan of handwritten notes). Can understand more intuitively by noting
		//  that the equation is:
		//
		// F1 = dEdv_parallel + 1/scalefactor(|v|) * [ cos( |v|/2 ) * dEdv_perpendicular + sin( |v|/2 ) * ( dEdv x u ) ]
		//
		// where u is the unit vector along v, v/|v|, and
		// where scalefactor( |v| ) =  sin( |v|/2 ) / (|v|/2).
		//
		// This is the same as the equation that rotates an arbitrary vector a by an angle theta around axis u:
		//
		// a_rotated = a_parallel + cos( theta ) * a_perpendicular + sin( theta ) * u x a
		//
		// where a_parallel = dot(a,u) u  [component of a parallel to unit vector u] and
		//       a_perpendicular = a - a_parallel [ component of a perpendicular to u] and
		//       a x u = is the other component of a perpendicular to both of the above.
		//
		//  [see, e.g., Staley, Quaternions and the Dirac Belt Trick, https://arxiv.org/abs/1001.1778 ]
		//
		// So the expression for F1 involves scaling down each component of dEdv that is perpendicular to v by
		//    sin( |v|/2 ) / (|v|/2),  and then reversing the rotation |v|/2, i.e. the angle theta of the quaternion
		//
		// The scalefactor is not that mysterious -- we know that the phase space volume in rotation vector representation
		//   is not constant but involves scaling, e.g. if we were integrating we'd use volume element
		//
		//   [sin( |v|/2) / (|v|/2)]^2 dv_x dv_y dv_z = dv_parallel  ( sin( |v|/2 )/ |v|/2) dv_perp1  ( sin( |v|/2 )/ |v|/2) dv_perp2
		//
		// Of course the whole thing is hard to derive completely intuitively unless you happen to be
		//  deeply familar with SO(3) and its embedding in R^4 via quaternions (see that Staley article).
		//
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////

		Vector const V_norm = rotvector / V;
		f1_rot += dot( dEdv, V_norm ) * V_norm;
		// note: cot( x ) = tan( pi/2 - x )
		f1_rot += (V/2) * std::tan( pi/2 - V/2 )*( dEdv - dot( dEdv , V_norm )*V_norm );
		f1_rot += -(V/2) * cross( dEdv, V_norm);
	} else {
		f1_rot += dEdv;
	}
	// Rotating into global frame from local coordinate frame at takeoff atom -- same
	//  procedure as above regarding translational contribution.
	Vector const & f1_rot_global = stub1_.rotation() * f1_rot;
	f1 += f1_rot_global;

	/////////////////////////////////////////////////////////////////////
	// For the landing atom, f1, f2 are negative of above.
	// Note that the f1/f2 have to cancel exactly
	//  from landing atom and takeoff atom, when accumulation of F1 and F2
	//  down the atom tree passes through, e.g., the takeoff atom and then
	//  later through the landing atom, to other atoms 'outside' these
	//  two atoms whose torsions cannot impact the energy via this
	//  potential.
	// Easy to check the similarity of downstream vs. upstream explicitly:
	// https://github.com/rhiju/loop_close/blob/master/notes/Das_6DPotential_F1_F2.pdf
	/////////////////////////////////////////////////////////////////////
	if ( !takeoff ) {
		f1 *= -1.0;
		f2 *= -1.0;
	}

}

} //evaluator
} //loop_graph
} //scoring
} //core
