// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/Jump.cc
/// @brief  Kinematics Jump class
/// @author Phil Bradley


// Unit headers
#include <core/kinematics/Jump.hh>

// Package headers
#include <core/kinematics/Stub.hh>

// Numeric Headers
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.hh>

// C++ Headers
#include <iostream>

#include <utility/vector1.hh>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>

#endif // SERIALIZATION


namespace core {
namespace kinematics {


//static const utility::vector1<double> ZERO( 6, 0.0 );

// convenience
numeric::xyzMatrix< core::Real >
Jump::mirror_z_transform = numeric::xyzMatrix< core::Real >::rows( 1.,0.,0., 0.,1.,0., 0.,0.,-1. );


///////////////////////////////////////////////////////////////////////////////
// reset
void
Jump::reset()
{
	rt_.reset();
	rb_delta[1] = rb_delta[2] = ZERO;
	rb_center[1] = rb_center[2] = Vector(0.0);
	invert_downstream_ = invert_upstream_ = false;
}


///////////////////////////////////////////////////////////////////////////////
bool
Jump::ortho_check() const
{
	return rt_.ortho_check();
}


///////////////////////////////////////////////////////////////////////////////
void
Jump::fold_in_rb_deltas()
{
	// n2c
	rt_.fold_in_rb_deltas( rb_delta[1], rb_center[1] );
	rt_.reverse();

	// c2n
	rt_.fold_in_rb_deltas( rb_delta[2], rb_center[2] );
	rt_.reverse();

	rb_delta[1] = rb_delta[2] = ZERO;
}

void
Jump::set_invert( bool upstream, bool downstream ) {
	invert_upstream_ = upstream;
	invert_downstream_ = downstream;
}

///////////////////////////////////////////////////////////////////////////////
/// @details stub defines the local reference frame for the corresponding direction
/// of the jump (forward == folding direction, backward = reverse)
/// by which to transform center\n
/// center is an xyz position in the absolute (protein) reference frame\n
/// the column vectors of stub are unit vectors defined in the
/// absolute reference frame( the frame in which center is defined)\n
/// so left-multiplying by stub.M.transposed() puts an absolute ref
/// frame vector (like center) into the local jump reference frame\n
/// dir is the direction through the jump,
///  1 == forward == folding order,
/// -1 == backward == reverse folding order\n
/// if dir == 1, the stub should be the stub for the downstream
/// jump atom, (eg returned by "folder.downstream_jump_stub( jump_number )" or
/// "downstream_jump_atom.get_stub()" )
/// and the center is the center of rotation for the
/// downstream segment (ie the segment which is built later during folding).
/// The center determines how the forward rb_deltas
/// are interpreted, and is the center for rigid-body rotation during
/// minimization.\n
/// if dir == -1, the stub should be the stub for the upstream
/// jump atom, (eg returned by folder.upstream_jump_stub( jump_number ) )
/// and the center is the center of rotation for the
/// upstream segment. The center determines how the reverse rb_deltas
/// are interpreted, eg in perturbation moves when we modify the
/// reverse rb-deltas directly, eg in the call
/// "Jump::gaussian_move(-1,tmag,rmag)"\n
///
void
Jump::set_rb_center(
	const int dir,
	Stub const & stub,
	Vector const & center
)
{
	fold_in_rb_deltas();
	int const index( dir == 1 ? 1 : 2 );
	rb_center[index] = stub.M.transposed() * ( center - stub.v );
}
///////////////////////////////////////////////////////////////////////////////
bool
Jump::nonzero_deltas() const
{
	const Real tol(1.0e-3);
	Real sum(0.0);
	for ( int i=1; i<=2; ++i ) {
		for ( int j=1; j<= 6; ++j ) {
			sum += std::abs( rb_delta[i][j] );
		}
	}
	return ( sum > tol );
}


///////////////////////////////////////////////////////////////////////////////
/// @details choose a random unit vector as the direction of translation,
/// by selecting its spherical coordinates phi(0-180) and theta(0-360)
/// then move dist_in along that direction
void
Jump::random_trans( const float dist_in )
{

	using numeric::y_rotation_matrix_degrees;
	using numeric::z_rotation_matrix_degrees;
	using numeric::conversions::degrees;
	using numeric::sin_cos_range;

	const Real theta( 360.0 * numeric::random::rg().uniform());
	const Real phi( degrees( std::acos(sin_cos_range(1.0-2.0*numeric::random::rg().uniform()))));
	const Real dist( dist_in );

	fold_in_rb_deltas();
	rt_.set_translation( dist * (
		y_rotation_matrix_degrees(phi) *
		z_rotation_matrix_degrees(theta) ).col_z() );
}

///////////////////////////////////////////////////////////////////////////////
/// @details clear existing rb_delta first if any and then apply the gaussian move
/// Return the move that was applied.
utility::vector1<Real>
Jump::gaussian_move(int const dir, float const trans_mag, float const rot_mag) {
	fold_in_rb_deltas(); // clear rb_delta
	for ( int i = 1; i <= 3; ++i ) {
		set_rb_delta( i, dir, Real( trans_mag * numeric::random::rg().gaussian() ) );
		set_rb_delta( i+3, dir, Real( rot_mag * numeric::random::rg().gaussian() ) );
	}
	utility::vector1<Real> this_rb_delta =  get_rb_delta(dir);
	fold_in_rb_deltas();
	return this_rb_delta;
}

/// @details do a gaussian move along a certain rb direction
void
Jump::gaussian_move_single_rb(
	int const dir,
	float const mag,
	int rb
)
{
	fold_in_rb_deltas(); // clear rb_delta
	set_rb_delta( rb, dir, Real( mag * numeric::random::rg().gaussian() ) );
	fold_in_rb_deltas();
}

///////////////////////////////////////////////////////////////////////////////
// set rb_delta
void
Jump::set_rb_delta(
	int const rb_no,
	int const dir,
	Real const value
)
{
	rb_delta[ rb_index(dir) ][rb_no] = value;
}

/// @brief set rb_deltas by direction
void
Jump::set_rb_deltas( int const dir, utility::vector1<Real> const & rb_deltas){
	rb_delta[ rb_index(dir) ] = rb_deltas;
}


///////////////////////////////////////////////////////////////////////////////
void
Jump::set_rotation(
	Matrix const & R
)
{
	fold_in_rb_deltas(); // clear rb_delta
	rt_.set_rotation( R );
}


///////////////////////////////////////////////////////////////////////////////
void
Jump::set_translation(
	Vector const & t
)
{
	fold_in_rb_deltas(); // clear rb_delta
	rt_.set_translation( t );
}


///////////////////////////////////////////////////////////////////////////////
/// @details stub should be the Stub of the upstream jump.
/// This routine has the effect (I think) of rotating the downstream stub
/// by the rotation "matrix" about the center "center", both of which are written
/// in xyz lab frame.
void
Jump::rotation_by_matrix(
	Stub const & stub,
	Vector const & center, //in xyz (absolute- or protein-reference-) frame
	Matrix const & matrix //in xyz frame
)
{
	fold_in_rb_deltas(); // clear rb_delta

	Matrix const & m( stub.M );

	// find the rotation center in jump coordinate system
	Vector new_center = m.transposed() * ( center - stub.v);

	// find the operator matrix when transofrmed to jump coord system
	Matrix new_matrix = m.transposed() * ( matrix * m );

	// new translation vector after applying matrix in jump coord system
	rt_.set_translation( new_center + new_matrix * ( rt_.get_translation() -
		new_center ) );

	// new rotation after applying matrix in jump coord system
	rt_.set_rotation( new_matrix * rt_.get_rotation() );
}

////////////////////////////////////////////////////////////////////////////////
/// @details stub should be the Stub of the upstream jump.
/// this does a rotation by alpha degrees around the specified
/// the axis and center, both of which are written in xyz frame
void
Jump::rotation_by_axis(
	Stub const & stub,
	Vector const & axis,
	Vector const & center, //in xyz frame
	float const alpha ///< degrees
)
{
	using numeric::conversions::radians;

	fold_in_rb_deltas();
	Matrix mat = numeric::rotation_matrix( axis, Real( radians( alpha ) ) );
	rotation_by_matrix( stub, center, mat );
}


///////////////////////////////////////////////////////////////////////////////
/// @details stub should be the Stub of the upstream jump.
/// this does a translation along the axis by dist.
/// the axis is  written in xyz frame

void
Jump::translation_along_axis(
	Stub const & stub,
	Vector const & axis, //in xyz frame
	float const dist // in angstrom
)
{
	fold_in_rb_deltas();
	Vector new_trans( Real(dist) * axis.normalized() );
	rt_.set_translation( rt_.get_translation() + stub.M.transposed()*new_trans );
}


///////////////////////////////////////////////////////////////////////////////
///
void
Jump::reverse()
{
	fold_in_rb_deltas();
	rt_.reverse();
}

///////////////////////////////////////////////////////////////////////////////
///
Jump
Jump::reversed() const
{
	Jump ret_val( *this );
	ret_val.fold_in_rb_deltas();
	ret_val.reverse();
	return ret_val;
}

///////////////////////////////////////////////////////////////////////////////
void
Jump::identity_transform()
{
	rb_delta[1] = rb_delta[2] = ZERO;
	rb_center[1] = rb_center[2] = Vector(0.0);
	rt_.identity_transform();
}


///////////////////////////////////////////////////////////////////////////////
/// @details fold all the rb_delta first and then make the jump.
/// make jump with stubs instead of FArrays
void
Jump::make_jump(
	Stub const & stub1,
	Stub & stub2
) const
{
	// make temporary local copy of our rotation-translation
	RT tmp_rt(rt_);

	// n2c
	tmp_rt.fold_in_rb_deltas( rb_delta[1], rb_center[1] );
	tmp_rt.reverse();

	// c2n
	tmp_rt.fold_in_rb_deltas( rb_delta[2], rb_center[2]  );

	// back to original direction
	tmp_rt.reverse();

	// main logic for inversion lives here
	// ensure correct handedness of input stub
	//bool stub1_inverted = (stub1.M.det() < 0);
	//runtime_assert( stub1_inverted == invert_upstream_ ); //fd during pose construction this is not true
	tmp_rt.make_jump( stub1, stub2 );

	// is this jump an inverting jump?
	//bool stub2_inverted = (stub2.M.det() < 0);

	if ( invert_upstream_ != invert_downstream_ ) {
		stub2.M = stub2.M * Jump::mirror_z_transform;
	}
}


///////////////////////////////////////////////////////////////////////////////
/// @note: we dont reset rb_center!!!!!!!!!
void
Jump::from_stubs(
	Stub const & stub1,
	Stub const & stub2
)
{
	// here we ignore inversion.  R/T of a jump is always stored w.r.t. a right-handed system
	rt_.from_stubs( stub1, stub2 );
	rb_delta[1] = rb_delta[2] = ZERO;
}

///////////////////////////////////////////////////////////////////////////////
/// @details this function translate a virtual bond type constraints into jump
///transformation.
///in atoms vector, A1, A2, A3, B1, B2, B3,
///and A1-B1 forms the virtual bond, like A3-A2-A1---B1-B2-B3
///in csts vector, disAB, angleA and B, didedral A, AB and B.
///B1, B2 and B3 are moved given cst definition.
void
Jump::from_bond_cst(
	utility::vector1< Vector > & atoms,
	utility::vector1< Real > const & csts
)
{
	debug_assert( atoms.size() == 6 && csts.size() == 6 );
	// make sure A3-A2-A1 or B3-B2-B1 is not linear (i.e., a stub can be defined)
	Stub stubA(atoms[1], atoms[2], atoms[3]); // a1, a2, a3
	Stub stubB(atoms[4], atoms[5], atoms[6]); // b1, b2, b3
	debug_assert( stubA.is_orthogonal(1e-3) && stubB.is_orthogonal(1e-3) );
	// udpate B1, B2, B3 positions
	Vector b1 = stubA.spherical( csts[4]/*dihedralA*/, csts[2]/*angleA*/, csts[1]/*disAB*/ );
	Stub stub_b1( b1, atoms[1], atoms[2]); // b1, a1, a2
	debug_assert(stub_b1.is_orthogonal(1e-3));
	Vector b2 = stub_b1.spherical( csts[5]/*dihedralAB*/, csts[3]/*angleB*/, atoms[4].distance(atoms[5]) );
	Stub stub_b2( b2, b1, atoms[1]); // b2, b1, a1
	debug_assert(stub_b2.is_orthogonal(1e-3));
	Vector b3 = stub_b2.spherical( csts[6]/*dihedralB*/, angle_of(atoms[4], atoms[5], atoms[6]), atoms[5].distance( atoms[6] ) );
	// update new b1, b2, b3 and get new  RT
	Stub new_stubB(b1,b2,b3);
	debug_assert( new_stubB.is_orthogonal(1e-3) );
	rt_.from_stubs( stubA, new_stubB );
	rb_delta[1] = rb_delta[2] = ZERO;
	atoms[4] = b1;
	atoms[5] = b2;
	atoms[6] = b3;
	return;
}
///////////////////////////////////////////////////////////////////////////////
// input
std::istream &
operator >>(
	std::istream & is,
	Jump & jump
)
{
	is >> jump.rt_;
	jump.rb_delta[1] = jump.rb_delta[2] = ZERO;
	jump.rb_center[1] = jump.rb_center[2] = Jump::Vector( 0.0 );
	return is;
}


///////////////////////////////////////////////////////////////////////////////
// stream output:
std::ostream&
operator <<(
	std::ostream & os,
	const Jump & jump
)
{
	if ( jump.get_invert_upstream() || jump.get_invert_downstream() ) {
		os << "mirrored: " << jump.get_invert_upstream() << "/" << jump.get_invert_downstream() << " ";
	}
	if ( jump.nonzero_deltas() ) {
		// old-style verbose output
		os << "Jump:: nonzero_deltas= " << jump.nonzero_deltas() << '\n';
		os << jump.rt_;
		for ( int i = 1; i <= 2; ++i ) {
			std::string tag ( (i==1) ? "n2c" : "c2n" );
			os << "rb_delta " << tag;
			for ( int j = 1; j <= 6; ++j ) {
				//os << F(9,3,jump.rb_delta(j,i));
				os << jump.rb_delta[i][j] << ' '; //chu -- more digits for precision
			}
			os << '\n';
			os << "rb_center " << tag;
			for ( int j = 0; j <= 2; ++j ) { //fpd Vector is 0-indexed
				os << jump.rb_center[i][j] << ' ';
			}
			os << '\n';
		}
	} else {
		//
		os << jump.rt_;
	}
	return os;
}

///////////////////////////////////////////////////////////////////////////////
Real
distance(
	Jump const & a_in,
	Jump const & b_in
)
{

	Jump a( a_in ), b( b_in );

	a.fold_in_rb_deltas();
	b.fold_in_rb_deltas();

	return distance( a.rt_, b.rt_ );
}


///////////////////////////////////////////////////////////////////////////////
void
jump_distance(
	Jump const & a_in,
	Jump const & b_in,
	Real & dist,
	Real & theta
)
{

	Jump a( a_in ), b( b_in );

	a.fold_in_rb_deltas();
	b.fold_in_rb_deltas();

	Jump::Matrix A( a.get_rotation() );
	Jump::Vector v( a.get_translation() );

	Jump::Matrix B( b.get_rotation() );
	Jump::Vector w( b.get_translation() );

	//Real theta; // in radians
	// A * B^-1 = rotation to make A==B.
	rotation_axis( A * B.transposed(), theta );

	dist = v.distance(w);

	//return distance( v, w ) + theta;
}

RT Jump::rt() const {
	if ( !nonzero_deltas() ) {
		return rt_;
	} else {
		Jump tmp( *this );
		tmp.fold_in_rb_deltas();
		return tmp.rt_;
	}
}

bool Jump::operator==( Jump const& other ) const {
	if ( this == &other ) {
		return true;
	} else if ( this->get_rotation() == other.get_rotation() &&
			this->get_translation() == other.get_translation() ) {
		return true;
	} else {
		return false;
	}
}

} // namespace kinematics
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::kinematics::Jump::save( Archive & arc ) const {
	arc( CEREAL_NVP( rt_ ) ); // class core::kinematics::RT
	arc( CEREAL_NVP( rb_delta ) ); // utility::vector1<utility::vector1<Real> >
	arc( CEREAL_NVP( rb_center ) ); // utility::vector1<Vector>
	arc( CEREAL_NVP( invert_upstream_ ) ); // _Bool
	arc( CEREAL_NVP( invert_downstream_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::kinematics::Jump::load( Archive & arc ) {
	arc( rt_ ); // class core::kinematics::RT
	arc( rb_delta ); // utility::vector1<utility::vector1<Real> >
	arc( rb_center ); // utility::vector1<Vector>
	arc( invert_upstream_ ); // _Bool
	arc( invert_downstream_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( core::kinematics::Jump );
#endif // SERIALIZATION
