// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/xyzTransform.hh
/// @brief  Fast rigid xform 3x3 matrix + vector
/// @author will sheffler


#ifndef INCLUDED_numeric_xyzTransform_hh
#define INCLUDED_numeric_xyzTransform_hh


// Unit headers
#include <numeric/xyzTransform.fwd.hh>

// Package headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>
#include <utility/vector1.hh>
// #include <numeric/xyzVector.io.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/numbers.hh>


// C++ headers
#include <utility/assert.hh>
#include <cmath>


namespace numeric {

template< typename T >
class xyzTransform {
public:
	typedef xyzTransform<T>  Transform;
	typedef xyzMatrix<T>     Matrix;
	typedef xyzVector<T>     Vector;
	typedef xyzVector<T>     V;

	typedef utility::fixedsizearray1< T, 6 > T6;

	typedef std::numeric_limits<T> LIM;

	Matrix R;
	Vector t;

	xyzTransform() : R(xyzMatrix<T>::identity()),t(0,0,0) {}
	xyzTransform(Matrix const & rin) : R(rin),t(0,0,0) {}
	xyzTransform(Vector const & tin) : R(xyzMatrix<T>::identity()),t(tin) {}
	xyzTransform(Matrix const & rin, Vector const & tin) : R(rin),t(tin) {}
	xyzTransform(T6 const & _rt6){ rt6(_rt6); }


	Transform & from_four_points(Vector const & c, Vector const & u, Vector const & v, Vector const & w){
		Vector e1( u - v);
		e1.normalize();
		Vector e3( cross( e1, w - v ) );
		e3.normalize();
		Vector e2( cross( e3,e1) );
		R.col_x( e1 ).col_y( e2 ).col_z( e3 );
		t = c;
		return *this;
	}
	xyzTransform(Vector const & u, Vector const & v, Vector const & w) { from_four_points(u,u,v,w); }
	xyzTransform(Vector const & c, Vector const & u, Vector const & v, Vector const & w) { from_four_points(c,u,v,w); }

	static xyzTransform<T> identity() { return xyzTransform(); }
	static xyzTransform<T> BAD_XFORM() { return xyzTransform( Matrix::cols(LIM::quiet_NaN(),LIM::quiet_NaN(),LIM::quiet_NaN(),LIM::quiet_NaN(),LIM::quiet_NaN(),LIM::quiet_NaN(),LIM::quiet_NaN(),LIM::quiet_NaN(),LIM::quiet_NaN()), Vector(LIM::quiet_NaN(),LIM::quiet_NaN(),LIM::quiet_NaN()) ); }
	static T6 BAD_RT6() { return T6(LIM::quiet_NaN()); }

	T const & xx() const { return R.xx(); }
	T const & xy() const { return R.xy(); }
	T const & xz() const { return R.xz(); }
	T const & yx() const { return R.yx(); }
	T const & yy() const { return R.yy(); }
	T const & yz() const { return R.yz(); }
	T const & zx() const { return R.zx(); }
	T const & zy() const { return R.zy(); }
	T const & zz() const { return R.zz(); }
	T const & px() const { return t.x (); }
	T const & py() const { return t.y (); }
	T const & pz() const { return t.z (); }
	T const &  x() const { return t.x (); }
	T const &  y() const { return t.y (); }
	T const &  z() const { return t.z (); }
	T       & xx()       { return R.xx(); }
	T       & xy()       { return R.xy(); }
	T       & xz()       { return R.xz(); }
	T       & yx()       { return R.yx(); }
	T       & yy()       { return R.yy(); }
	T       & yz()       { return R.yz(); }
	T       & zx()       { return R.zx(); }
	T       & zy()       { return R.zy(); }
	T       & zz()       { return R.zz(); }
	T       & px()       { return t.x (); }
	T       & py()       { return t.y (); }
	T       & pz()       { return t.z (); }
	T       &  x()       { return t.x (); }
	T       &  y()       { return t.y (); }
	T       &  z()       { return t.z (); }


	Transform inverse          () const { return xyzTransform(R.transposed(),R.transposed()*-t); }
	Transform operator~       (       ) const { return xyzTransform(R.transposed(),R.transposed()*-t); }
	T distance        ( Transform const & b ) const { return std::sqrt( distance_squared(b) ); }
	T distance_squared( Transform const & b ) const { return R.col_x().distance_squared(b.R.col_x() ) + R.col_y().distance_squared(b.R.col_y() ) +
		R.col_z().distance_squared(b.R.col_z() ) + t.distance_squared( b.t ); }
	friend bool operator==( Transform const & a, Transform const & b ){ return a.distance_squared(b)<=0.000001; }
	friend bool operator!=( Transform const & a, Transform const & b ){ return a.distance_squared(b)> 0.000001; }

	friend    Transform operator +( Transform const & a, Vector const & b ){ return Transform( a.R, a.t+a.R*b ); }
	friend    Transform operator -( Transform const & a, Vector const & b ){ return Transform( a.R, a.t-a.R*b ); }
	friend    Transform operator +( Vector const & a, Transform const & b ){ return Transform( b.R, b.t+    a ); }
	friend    Transform operator -( Vector const & a, Transform const & b ){ return Transform( b.R, b.t-    a ); }
	friend    Transform operator *( Transform const & a, T const & b ){ return Transform( a.R, a.t*b ); }
	friend    Transform operator *( T const & a, Transform const & b ){ return Transform( b.R, b.t*a ); }
	friend    Transform operator *( Transform const & a, Matrix const & b ){ return Transform( a.R*b, a.t ); }
	friend    Transform operator *( Matrix const & a, Transform const & b ){ return Transform( a*b.R, a*b.t ); }
	friend    Vector operator *( Transform const & x, Vector const & v ){ return x.R*v+x.t; }
	friend    Transform operator *( Transform const & a, Transform const & b ){ return Transform( a.R*b.R, a.R*b.t + a.t ); }

	// friend    Transform operator /( Transform const & n, Transform const & d ){ return (~d)*n; }
	friend    Transform operator -( Transform const & a, Transform const & b ){ return Transform( b.R.transposed()*a.R, b.R.transposed() * (a.t-b.t) ); }
	static    Transform rot       ( Matrix const & rot, Vector const & o_cen, Vector const & cen){ return xyzTransform<T>( rot, rot*-o_cen + cen);} //Correct Version, Translates to Origin based on original centroid
	static    Transform rot       ( Matrix const & rot, Vector const & cen = Vector(0,0,0)){ return Transform( rot, rot*-cen + cen );}
	static    Transform rot       ( Vector const & axs, T const & ang, Vector const & cen = Vector(0,0,0)){ return rot(rotation_matrix        (axs,ang),cen); }
	static    Transform rot_deg   ( Vector const & axs, T const & ang, Vector const & cen = Vector(0,0,0)){ return rot(rotation_matrix_degrees(axs,ang),cen); }

	static    Transform align( Vector const & to, Vector const & from = Vector(1,0,0) ){
		Vector  const axis = to.cross(from);
		T const angle = std::acos(to.normalized().dot(from.normalized()));
		return Transform(rotation_matrix(axis,-angle));
	}
	static    Transform align_fast( V const & to, V const & from = V(1.0, 0.0, 0.0) ){ // rot 180 around avg. axis
		return Transform( T(2.0)*projection_matrix<T>(to+from)-xyzMatrix<T>::identity() );
	}

	Vector     xform( Vector const & v ) const { return R*v+t; }
	Vector inv_xform( Vector const & v ) const { return R.transposed()*(v-t); }
	template<typename T2> T2 operator()(T2 const & x) { return (*this)*x; }

	void to_quaternion( T& qw, T& qx, T& qy, T& qz ) const {
		T const t = xx()+yy()+zz(); // trace of M
		T const r = sqrt(1+t);
		qw = 0.5*r;
		qx = copysign(0.5*sqrt(1+xx()-yy()-zz()), yz()-zy());
		qy = copysign(0.5*sqrt(1-xx()+yy()-zz()), zx()-xz());
		qz = copysign(0.5*sqrt(1-xx()-yy()+zz()), xy()-yx());
		qz = qx>0 ? qz : -qz;
		qy = qx>0 ? qy : -qy;
		qw = qx>0 ? qw : -qw;
		qx = qx>0 ? qx : -qx; // qx must be last!
		// double const t = xx()+yy()+zz();
		// double const r = sqrt(1+t);
		// double const s = 0.5/r;
		// double const qw = 0.5*r;
		// double const qx = (yz()-zy())*s;
		// double const qy = (zx()-xz())*s;
		// double const qz = (xy()-yx())*s;
	}
	void from_quaternion( T const & qw, T const & qx, T const & qy, T const & qz ){
		// X x( numeric::xyzMatrix<double>::cols(
		//     1.0 - 2.0*qy*qy - 2.0*qz*qz,       2.0*qx*qy - 2.0*qz*qw,       2.0*qx*qz + 2.0*qy*qw,
		//           2.0*qx*qy + 2.0*qz*qw, 1.0 - 2.0*qx*qx - 2.0*qz*qz,       2.0*qy*qz - 2.0*qx*qw,
		//           2.0*qx*qz - 2.0*qy*qw,       2.0*qy*qz + 2.0*qx*qw, 1.0 - 2.0*qx*qx - 2.0*qy*qy   ));
		T const Nq = qw*qw + qx*qx + qy*qy + qz*qz;
		T const s  = (Nq > 0.0) ? 2.0/Nq : 0.0;
		T const X  = qx*s; T const Y  = qy*s; T const Z  = qz*s;
		T const wX = qw*X; T const wY = qw*Y; T const wZ = qw*Z;
		T const xX = qx*X; T const xY = qx*Y; T const xZ = qx*Z;
		T const yY = qy*Y; T const yZ = qy*Z; T const zZ = qz*Z;
		xx() = 1.0-(yY+zZ); yx() = xY-wZ;       zx() = xZ+wY;
		xy() = xY+wZ;  yy() = 1.0-(xX+zZ); zy() = yZ-wX;
		xz() = xZ-wY;  yz() = yZ+wX;        zz() = 1.0-(xX+yY);
	}

	/// @brief see numeric/HomogeneousTransform
	// guarantee 0 <= e(1) < pi && 0 <= e(2) < pi && 0 <= e(3) <= pi
	xyzVector< T > euler_angles_rad() const {
		xyzVector< T > euler;
		T const FLOAT_PRECISION( 1e-5 );
		if ( R.zz() >= 1 - FLOAT_PRECISION ) {
			euler(1) = std::atan2( sin_cos_range( R.yx() ), sin_cos_range( R.xx() ) );
			euler(2) = 0.0;
			euler(3) = 0.0;
		} else if ( R.zz() <= -1 + FLOAT_PRECISION ) {
			euler(1) = std::atan2( sin_cos_range( R.yx() ), sin_cos_range( R.xx() ) );
			euler(2) = 0.0;
			euler(3) = (T)constants::d::pi;
		} else {
			T pos_sin_theta = std::sqrt( 1 - R.zz()*R.zz() ); // sin2theta = 1 - cos2theta.
			euler(3) = std::asin( pos_sin_theta );
			if ( R.zz() < 0 ) {
				euler(3) = (T)constants::d::pi - euler(3);
			}
			euler(1) = std::atan2( R.xz(), -R.yz() );
			euler(2) = std::atan2( R.zx(),  R.zy() );
		}
		euler(1) += euler(1)<0.0 ? (T)constants::d::pi_2 : 0.0;
		euler(2) += euler(2)<0.0 ? (T)constants::d::pi_2 : 0.0;

		euler(1) = min<T>(max<T>(0.0,euler(1)),(T)constants::d::pi_2-0.000000000001);
		euler(2) = min<T>(max<T>(0.0,euler(2)),(T)constants::d::pi_2-0.000000000001);
		euler(3) = min<T>(max<T>(0.0,euler(3)),(T)constants::d::pi  -0.000000000001);

		assert( 0 <= euler(1) ); assert( euler(1) <  (T)constants::d::pi_2 );
		assert( 0 <= euler(2) ); assert( euler(2) <  (T)constants::d::pi_2 );
		assert( 0 <= euler(3) ); assert( euler(3) <= (T)constants::d::pi   );
		return euler;
	}

	// xyzVector< T > euler_angles_rad_fast() const {
	//  xyzVector< T > euler;
	//  T const FLOAT_PRECISION( 1e-5 );
	//  if ( R.zz() >= 1 - FLOAT_PRECISION ){
	//   euler(1) = fastatan2( sin_cos_range( R.yx() ), sin_cos_range( R.xx() ) );
	//   euler(2) = 0.0;
	//   euler(3) = 0.0;
	//  } else if ( R.zz() <= -1 + FLOAT_PRECISION ){
	//   euler(1) = fastatan2( sin_cos_range( R.yx() ), sin_cos_range( R.xx() ) );
	//   euler(2) = 0.0;
	//   euler(3) = (T)constants::d::pi;
	//  } else {
	//   T pos_sin_theta = std::sqrt( 1 - R.zz()*R.zz() ); // sin2theta = 1 - cos2theta.
	//   euler(3) = fastasin( pos_sin_theta );
	//   if ( R.zz() < 0 ) {
	//    euler(3) = (T)constants::d::pi - euler(3);
	//   }
	//   euler(1) = fastatan2( R.xz(), -R.yz() );
	//   euler(2) = fastatan2( R.zx(),  R.zy() );
	//  }
	//  euler(1) += euler(1)<0.0 ? (T)constants::d::pi_2 : 0.0;
	//  euler(2) += euler(2)<0.0 ? (T)constants::d::pi_2 : 0.0;
	//  euler(1) = min(max(0.0,euler(1)),(T)constants::d::pi_2-0.000000000001);
	//  euler(2) = min(max(0.0,euler(2)),(T)constants::d::pi_2-0.000000000001);
	//  euler(3) = min(max(0.0,euler(3)),(T)constants::d::pi  -0.000000000001);
	//  assert( 0 <= euler(1) ); assert( euler(1) <  (T)constants::d::pi_2 );
	//  assert( 0 <= euler(2) ); assert( euler(2) <  (T)constants::d::pi_2 );
	//  assert( 0 <= euler(3) ); assert( euler(3) <= (T)constants::d::pi   );
	//  return euler;
	// }

	xyzVector< T > euler_angles_deg() const {
		return T(constants::d::radians_to_degrees) * (Vector)euler_angles_rad();
	}
	// xyzVector< T > euler_angles_deg_fast() const {
	//  return (T)constants::d::radians_to_degrees * (V)euler_angles_rad_fast();
	// }

	xyzTransform<T> & from_euler_angles_rad( T const & phi, T const & psi, T const & theta ) {
		// TC phi( euler( 1 ) ), psi( euler( 2 ) ), theta( euler( 3 ) );
		T const cos_phi( std::cos( phi ) ),    sin_phi( std::sin( phi ) );
		T const cos_psi( std::cos( psi ) ),    sin_psi( std::sin( psi ) );
		T const cos_theta( std::cos( theta )), sin_theta( std::sin( theta ) );
		R.xx() =  cos_psi * cos_phi - cos_theta * sin_phi * sin_psi;  R.yx() =  cos_psi * sin_phi + cos_theta * cos_phi * sin_psi;  R.zx() =  sin_psi * sin_theta;
		R.xy() = -sin_psi * cos_phi - cos_theta * sin_phi * cos_psi;  R.yy() = -sin_psi * sin_phi + cos_theta * cos_phi * cos_psi;  R.zy() =  cos_psi * sin_theta;
		R.xz() =                sin_theta * sin_phi;                  R.yz() =                     -sin_theta * cos_phi;            R.zz() =        cos_theta;
		return *this;
	}
	xyzTransform<T> & from_euler_angles_rad( xyzVector< T > const & euler ) {
		return from_euler_angles_rad(euler.x(),euler.y(),euler.z());
	}
	xyzTransform<T> & from_euler_angles_deg( T const & phi, T const & psi, T const & theta ) {
		T const & C = (T)constants::d::degrees_to_radians;
		return from_euler_angles_rad( C*phi, C*psi, C*theta );
	}
	xyzTransform<T> & from_euler_angles_deg( xyzVector< T > const & euler ) {
		xyzVector< T > euler_rad( euler );
		euler_rad *= (T)constants::d::degrees_to_radians;
		return from_euler_angles_rad( euler_rad );
	}

	T6 rt6() const {
		if ( bad() ) return BAD_RT6();
		T6 rt6;
		Vector const e( euler_angles_deg() );
		rt6[1] = px();
		rt6[2] = py();
		rt6[3] = pz();
		rt6[4] = e.x();
		rt6[5] = e.y();
		rt6[6] = e.z();
		return rt6;
	}

	// Real6 rt6fast() const {
	//  if(bad()) return BAD_RT6();
	//  Real6 rt6;
	//  V const e( euler_angles_deg_fast() );
	//  rt6[1] = px();
	//  rt6[2] = py();
	//  rt6[3] = pz();
	//  rt6[4] = e.x();
	//  rt6[5] = e.y();
	//  rt6[6] = e.z();
	//  return rt6;
	// }

	xyzTransform & rt6(T6 const & rt6){
		if ( utility::isnan(rt6[1]) ) { *this = BAD_XFORM(); return *this; }
		px() = rt6[1];
		py() = rt6[2];
		pz() = rt6[3];
		from_euler_angles_deg(rt6[4],rt6[5],rt6[6]);
		return *this;
	}

	xyzTransform & rt6(T const & i, T const & j, T const & k, T const & l, T const & m, T const & n){
		if ( utility::isnan(i) ) { *this = BAD_XFORM(); return *this; }
		px() = i;
		py() = j;
		pz() = k;
		from_euler_angles_deg(l,m,n);
		return *this;
	}

	// could be collisions!
	uint64_t hash64( T const & cw=0.1, T const & aw=360./1024. ) const {
		uint64_t x = (uint64_t)fabs((t.x()+cw/2.0)/cw);
		uint64_t y = (uint64_t)fabs((t.y()+cw/2.0)/cw);
		uint64_t z = (uint64_t)fabs((t.z()+cw/2.0)/cw);
		Vector e( euler_angles_deg() );
		uint64_t a = (uint64_t)fmod(e.x()+aw/2.0,360.0);
		uint64_t b = (uint64_t)fmod(e.y()+aw/2.0,360.0);
		uint64_t c = (uint64_t)fmod(e.z()+aw/2.0,360.0);
		assert(a < 1024);
		assert(b < 1024);
		assert(c < 1024);
		uint64_t k = x ^ y<<10 ^ z<<20 ^ a<<30 ^ b<<40 ^ c<<50;
		k ^= ((uint64_t)(x<0.0))<<60 ^ ((uint64_t)(y<0.0))<<61 ^ ((uint64_t)(z<0.0))<<62 ;
		return k;
	}

	uint64_t symhash64( T const & cw=0.1, T const & aw=360./1024. ) const {
		uint64_t h1 = hash64(cw,aw);
		uint64_t h2 = inverse().hash64(cw,aw);
		return h1 > h2 ? h1 : h2;
	}


	// this code should go elsewhere....
	// Copyright 2001 softSurfer, 2012 Dan Sunday
	// This code may be freely used and modified for any purpose
	// providing that this copyright notice is included with it.
	// SoftSurfer makes no warranty for this code, and cannot be held
	// liable for any real or imagined damage resulting from its use.
	// Users of this code must verify correctness for their application.
	struct Line { Vector P0,P1; };
	struct Plane { Vector n,V0; };

	// intersect3D_2Planes(): find the 3D intersection of two planes
	//    Input:  two planes Pn1 and Pn2
	//    Output: *L = the intersection line (when it exists)
	//    Return: 0 = disjoint (no intersection)
	//            1 = the two  planes coincide
	//            2 =  intersection in the unique line *L
	int intersect3D_2Planes( Plane Pn1, Plane Pn2, Line* L ) const {
		Vector   u = Pn1.n .cross( Pn2.n );
		T ax = (u.x() >= 0 ? u.x() : -u.x());
		T ay = (u.y() >= 0 ? u.y() : -u.y());
		T az = (u.z() >= 0 ? u.z() : -u.z());

		// test if the two planes are parallel
		if ( (ax+ay+az) < 0.000000001 ) {        // Pn1 and Pn2 are near parallel
			// test if disjoint or coincide
			Vector   v = Pn2.V0 -  Pn1.V0;
			if ( dot(Pn1.n, v) == 0 ) {          // Pn2.V0 lies in Pn1
				return 1;                    // Pn1 and Pn2 coincide
			} else {
				return 0;                    // Pn1 and Pn2 are disjoint
			}
		}

		// Pn1 and Pn2 intersect in a line
		// first determine max abs coordinate of cross product
		int maxc;
		if ( ax > ay ) {
			if ( ax > az ) maxc = 1;
			else         maxc = 3;
		} else {
			if ( ay > az ) maxc = 2;
			else         maxc = 3;
		}

		// next, to get a point on the intersect line
		// zero the max coord, and solve for the other two
		Vector    iP;                // intersect point
		T    d1, d2;            // the constants in the 2 plane equations
		d1 = -dot(Pn1.n, Pn1.V0);  // note: could be pre-stored  with plane
		d2 = -dot(Pn2.n, Pn2.V0);  // ditto

		switch (maxc) {             // select max coordinate
		case 1 :                     // intersect with x=0
			iP.x() = 0;
			iP.y() = (d2*Pn1.n.z() - d1*Pn2.n.z()) /  u.x();
			iP.z() = (d1*Pn2.n.y() - d2*Pn1.n.y()) /  u.x();
			break;
		case 2 :                     // intersect with y=0
			iP.x() = (d1*Pn2.n.z() - d2*Pn1.n.z()) /  u.y();
			iP.y() = 0;
			iP.z() = (d2*Pn1.n.x() - d1*Pn2.n.x()) /  u.y();
			break;
		case 3 :                     // intersect with z=0
			iP.x() = (d2*Pn1.n.y() - d1*Pn2.n.y()) /  u.z();
			iP.y() = (d1*Pn2.n.x() - d2*Pn1.n.x()) /  u.z();
			iP.z() = 0;
		}
		L->P0 = iP;
		L->P1 = iP + u;
		return 2;
	}


	void rotation_axis(Vector & axis, Vector & cen, T& angle) const {
		axis = numeric::rotation_axis<T>(R,angle);
		Vector const p1((T)-32.09501046777237,(T)03.36227004372687,(T) 35.34672781477340); // random...
		Vector const p2((T) 21.15113978202345,(T)12.55664537217840,(T)-37.48294301885574); // random...
		Vector const q1((*this)*p1);
		Vector const q2((*this)*p2);
		Vector const n1 = (q1-p1).normalized();
		Vector const n2 = (q2-p2).normalized();
		Vector const c1 = (p1+q1)/T(2.0);
		Vector const c2 = (p2+q2)/T(2.0);
		Plane Pn1 = {n1,c1};
		Plane Pn2 = {n2,c2};
		Line Linter;
		int inter_case = intersect3D_2Planes(Pn1,Pn2,&Linter);
		switch(inter_case){
		case 0 :
			cen = Vector(std::numeric_limits<T>::max(),std::numeric_limits<T>::max(),std::numeric_limits<T>::max());
			break;
		case 1 :
			cen = Vector(std::numeric_limits<T>::min(),std::numeric_limits<T>::min(),std::numeric_limits<T>::min());
			break;
		case 2 :
			Vector Laxis = (Linter.P1-Linter.P0).normalized();
			if ( -0.9999 < Laxis.dot(axis) && Laxis.dot(axis) < 0.9999 ) {
				// std::cout << Laxis << std::endl;
				// std::cout << axis  << std::endl;
				// std::cout << angle << std::endl;
				utility_exit_with_message("bad axis");
			}
			cen = Linter.P0;
			break;
		}
	}

	Vector rotation_axis() const {
		T angle;
		return numeric::rotation_axis<T>(R,angle);
	}

	T rotation_angle_degrees() const {
		T angle;
		numeric::rotation_axis<T>(R,angle);
		return numeric::conversions::degrees(angle);
	}
	T rotation_angle() const {
		T angle;
		numeric::rotation_axis<T>(R,angle);
		return angle;
	}

	T rotation_cosine() const {
		return sin_cos_range( ( R.trace() - (T)1.0 ) / (T)2.0 );
	}

	T rotation_sine() const {
		T cos = sin_cos_range( ( R.trace() - (T)1.0 ) / (T)2.0 );
		return sqrt(1.0-cos*cos);
	}

	T approx_lever_distance( Transform const & o, T const & lever=1.0 ) const {
		Transform tmp = o * this->inverse();
		T ang = tmp.rotation_angle();
		return sqrt(t.distance_squared(o.t) + ang*ang*lever*lever);
	}

	bool bad() const {
		return utility::isnan(px()) || utility::isnan(py()) || utility::isnan(pz()) ||
			utility::isnan(xx()) || utility::isnan(xy()) || utility::isnan(xz()) ||
			utility::isnan(yx()) || utility::isnan(yy()) || utility::isnan(yz()) ||
			utility::isnan(zx()) || utility::isnan(zy()) || utility::isnan(zz()) ;
	}
	bool badfast() const {
		return utility::isnan(px()) || utility::isnan(xx());
	}

};

struct XformHash32 : std::unary_function<xyzTransform<float>,uint64_t> {
	uint64_t operator()(xyzTransform<float> const & xform) const {
		uint64_t const *x = (uint64_t const *)(&xform);
		return x[0]^x[1]^x[2]^x[3]^x[4]^x[5];
	}
};
struct XformHash64 : std::unary_function<xyzTransform<double>,uint64_t> {
	uint64_t operator()(xyzTransform<double> const & xform) const {
		uint64_t const *x = (uint64_t const *)(&xform);
		return x[0]^x[1]^x[2]^x[3]^x[4]^x[5]^x[6]^x[7]^x[8]^x[9]^x[10]^x[11];
	}
};


struct Xforms : public utility::vector1<xyzTransform<numeric::Real> > {
	Xforms()
	: utility::vector1<xyzTransform<numeric::Real> >() {}
	Xforms(xyzTransform<numeric::Real> x)
	: utility::vector1<xyzTransform<numeric::Real> >(1,x) {}
	Xforms(Size const & N, xyzTransform<numeric::Real> x = xyzTransform<numeric::Real>())
	: utility::vector1<xyzTransform<numeric::Real> >(N,x) {}
};


template<typename T, class OutputIterator>
void
expand_xforms(
	OutputIterator container,
	xyzTransform<T> const & G1,
	xyzTransform<T> const & G2,
	xyzTransform<T> const & /*G3*/,
	int N=5,
	Real r=9e9,
	xyzVector<T> const & test_point=xyzVector<T>(Real(1.0),Real(3.0),Real(10.0))
){
	xyzTransform<T> I;
	utility::vector1<xyzVector<T> > seenit;
	for ( int i0a = 0; i0a < 2; ++i0a ) { xyzTransform<T> x0a( (i0a==1 ? G1 : I) *   I );
		for ( int i0b = 0; i0b < 2; ++i0b ) { xyzTransform<T> x0b( (i0b==1 ? G2 : I) * x0a );
			for ( int i1a = 0; i1a < 2; ++i1a ) { xyzTransform<T> x1a( (i1a==1 ? G1 : I) * x0b );
				for ( int i1b = 0; i1b < 2; ++i1b ) { xyzTransform<T> x1b( (i1b==1 ? G2 : I) * x1a );
					for ( int i2a = 0; i2a < 2; ++i2a ) { xyzTransform<T> x2a( (i2a==1 ? G1 : I) * x1b );
						for ( int i2b = 0; i2b < 2; ++i2b ) { xyzTransform<T> x2b( (i2b==1 ? G2 : I) * x2a );
							for ( int i3a = 0; i3a < 2; ++i3a ) { xyzTransform<T> x3a( (i3a==1 ? G1 : I) * x2b );
								for ( int i3b = 0; i3b < 2; ++i3b ) { xyzTransform<T> x3b( (i3b==1 ? G2 : I) * x3a );
									for ( int i4a = 0; i4a < 2; ++i4a ) { xyzTransform<T> x4a( (i4a==1 ? G1 : I) * x3b );
										for ( int i4b = 0; i4b < 2; ++i4b ) { xyzTransform<T> x4b( (i4b==1 ? G2 : I) * x4a );
											for ( int i5a = 0; i5a < 2; ++i5a ) { xyzTransform<T> x5a( (i5a==1 ? G1 : I) * x4b );
												for ( int i5b = 0; i5b < 2; ++i5b ) { xyzTransform<T> x5b( (i5b==1 ? G2 : I) * x5a );
													for ( int i6a = 0; i6a < 2; ++i6a ) { xyzTransform<T> x6a( (i6a==1 ? G1 : I) * x5b );
														for ( int i6b = 0; i6b < 2; ++i6b ) { xyzTransform<T> x6b( (i6b==1 ? G2 : I) * x6a );
															for ( int i7a = 0; i7a < 2; ++i7a ) { xyzTransform<T> x7a( (i7a==1 ? G1 : I) * x6b );
																for ( int i7b = 0; i7b < 2; ++i7b ) { xyzTransform<T> x7b( (i7b==1 ? G2 : I) * x7a );
																	for ( int i8a = 0; i8a < 2; ++i8a ) { xyzTransform<T> x8a( (i8a==1 ? G1 : I) * x7b );
																		for ( int i8b = 0; i8b < 2; ++i8b ) { xyzTransform<T> x8b( (i8b==1 ? G2 : I) * x8a );
																			for ( int i9a = 0; i9a < 2; ++i9a ) { xyzTransform<T> x9a( (i9a==1 ? G1 : I) * x8b );
																				for ( int i9b = 0; i9b < 2; ++i9b ) { xyzTransform<T> x9b( (i9b==1 ? G2 : I) * x9a );
																					for ( int iAa = 0; iAa < 2; ++iAa ) { xyzTransform<T> xAa( (i7a==1 ? G1 : I) * x9b );
																						for ( int iAb = 0; iAb < 2; ++iAb ) { xyzTransform<T> xAb( (i7b==1 ? G2 : I) * xAa );
																							for ( int iBa = 0; iBa < 2; ++iBa ) { xyzTransform<T> xBa( (i8a==1 ? G1 : I) * xAb );
																								for ( int iBb = 0; iBb < 2; ++iBb ) { xyzTransform<T> xBb( (i8b==1 ? G2 : I) * xBa );
																									for ( int iCa = 0; iCa < 2; ++iCa ) { xyzTransform<T> xCa( (i9a==1 ? G1 : I) * xBb );
																										for ( int iCb = 0; iCb < 2; ++iCb ) { xyzTransform<T> xCb( (i9b==1 ? G2 : I) * xCa );
																											xyzTransform<T> const & x(xCb);
																											if ( x.t.length() > r ) continue;
																											xyzVector<T> test = x*test_point;
																											bool redundant(false);
																											for ( typename utility::vector1<xyzVector<T> >::const_iterator i = seenit.begin(); i != seenit.end(); ++i ) {
																												if ( i->distance_squared(test) < 0.001 ) redundant = true;
																											}
																											if ( redundant ) continue;
																											seenit.push_back(test);
																											*container++ = x;
																										}
																									}
																									if ( 1==N ) return;
																								}
																							}
																							if ( 2==N ) return;
																						}
																					}
																					if ( 3==N ) return;
																				}
																			}
																			if ( 4==N ) return;
																		}
																	}
																	if ( 5==N ) return;
																}
															}
															if ( 6==N ) return;
														}
													}
													if ( 7==N ) return;
												}
											}
											if ( 8==N ) return;
										}
									}
									if ( 9==N ) return;
								}
							}
							if ( 10==N ) return;
						}
					}
					if ( 11==N ) return;
				}
			}
			if ( 12==N ) return;
		}
	}
}

template<typename T, class OutputIterator>
void
expand_xforms(
	OutputIterator container,
	xyzTransform<T> const & G1,
	xyzTransform<T> const & G2,
	int N=5,
	Real r=9e9,
	xyzVector<T> const & test_point=xyzVector<T>(Real(1.0),Real(3.0),Real(10.0))
){
	expand_xforms(container,G1,G2,xyzTransform<T>(),N,r,test_point);
}

// PyRosetta concreate version
class Py_xyzTransform_double : public xyzTransform<double> {};

} // namespace numeric


#endif // INCLUDED_numeric_xyzTransform_HH
