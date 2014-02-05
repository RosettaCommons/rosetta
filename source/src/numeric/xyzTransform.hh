// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/xyzTransform.hh
/// @brief  Fast rigid xform 3x3 matrix + vector
/// @author will sheffler


#ifndef INCLUDED_numeric_xyzTransform_hh
#define INCLUDED_numeric_xyzTransform_hh


// Unit headers
#include <numeric/xyzTransform.fwd.hh>

// Package headers
#include <numeric/xyzTriple.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>
#include <utility/vector1.hh>

// C++ headers
#include <cassert>


namespace numeric {

template< typename T >
class xyzTransform {
public:
	typedef xyzTransform<T>  X;
	typedef xyzMatrix<T>     M;
	typedef xyzVector<T>     V;
	typedef X const   XC;
	typedef M const   MC;
	typedef V const   VC;
	typedef X const & XCR;
	typedef M const & MCR;
	typedef V const & VCR;
	typedef T const & TCR;

	M R;
	V t;

	xyzTransform() : R(xyzMatrix<T>::identity()),t(0,0,0) {}
	xyzTransform(MCR rin) : R(rin),t(0,0,0) {}
	xyzTransform(VCR tin) : R(xyzMatrix<T>::identity()),t(tin) {}
	xyzTransform(MCR rin, VCR tin) : R(rin),t(tin) {}
	xyzTransform(VCR u, VCR v, VCR w) { from_four_points(u,u,v,w); }
	xyzTransform(VCR c, VCR u, VCR v, VCR w) { from_four_points(c,u,v,w); }

	void from_four_points(VCR c, VCR u, VCR v, VCR w){
		V e1( u - v);
		e1.normalize();
		V e3( cross( e1, w - v ) );
		e3.normalize();
		V e2( cross( e3,e1) );
		R.col_x( e1 ).col_y( e2 ).col_z( e3 );
		t = c;
	}


	inline X operator~       (       ) const { return xyzTransform(R.transposed(),R.transposed()*-t); }
	inline T distance        ( XCR b ) const { return std::sqrt( distance_squared(b) ); }
	inline T distance_squared( XCR b ) const { return R.col_x().distance_squared(b.R.col_x() ) + R.col_y().distance_squared(b.R.col_y() ) +
	                                                  R.col_z().distance_squared(b.R.col_z() ) + t.distance_squared( b.t ); }
	friend inline bool operator==( XCR a, XCR b ){ return ( a.R==b.R && a.t==b.t ); }
	friend inline bool operator!=( XCR a, XCR b ){ return ( a.R!=b.R || a.t!=b.t );	}

	// friend inline    X operator +( XCR a, XCR b ){ return X( a.R+b.R, a.t+b.t ); }
	// friend inline    X operator -( XCR a, XCR b ){ return X( a.R-b.R, a.t-b.t ); }
	// friend inline    X operator *( XCR a, TCR b ){ return X( a.R*b, a.t*b); }
	// friend inline    X operator *( TCR a, XCR b ){ return X( a*b.R, a*b.t); }
	// friend inline    X operator /( XCR a, TCR b ){ return X( a.R/b, a.t/b); }

	friend inline    X operator +( XCR a, VCR b ){ return X( a.R, a.t+a.R*b ); }
	friend inline    X operator -( XCR a, VCR b ){ return X( a.R, a.t-a.R*b ); }
	friend inline    X operator +( VCR a, XCR b ){ return X( b.R, b.t+    a ); }
	friend inline    X operator -( VCR a, XCR b ){ return X( b.R, b.t-    a ); }
	friend inline    X operator *( XCR a, MCR b ){ return X( a.R*b, a.t ); }
	friend inline    X operator *( MCR a, XCR b ){ return X( a*b.R, a*b.t ); }
	friend inline    V operator *( XCR x, VCR v ){ return x.R*v+x.t; }
	friend inline    X operator *( XCR a, XCR b ){ return X( a.R*b.R, a.R*b.t + a.t ); }
	// friend inline    X operator /( XCR n, XCR d ){ return (~d)*n; }
	friend inline    X operator -( XCR a, XCR b ){ return X( b.R.transposed()*a.R, b.R.transposed() * (a.t-b.t) ); }
	static inline    X rot       ( MCR rot,VCR o_cen,VCR cen){ return X( rot, rot*-o_cen + cen);} //Correct Version, Translates to Origin based on original centroid
	static inline    X rot       ( MCR rot,          VCR cen){ return X( rot, rot*-cen + cen );} //INCORRECT VERSION 
	static inline    X rot       ( VCR axs, TCR ang, VCR cen){ return rot(rotation_matrix        (axs,ang),cen); }
	static inline    X rot_deg   ( VCR axs, TCR ang, VCR cen){ return rot(rotation_matrix_degrees(axs,ang),cen); }

	inline V     xform( VCR v ) const { return R*v+t; }
	inline V inv_xform( VCR v ) const { return R.transposed()*(v-t); }
	template<typename T2> inline T2 operator()(T2 const & x) { return (*this)*x; }

	/// @brief see numeric/HomogeneousTransform
	xyzVector< T > euler_angles_rad() const {
		xyzVector< T > euler;
		T const FLOAT_PRECISION( 1e-9 );
		if ( R.zz() >= 1 - FLOAT_PRECISION ){
			euler(1) = std::acos( sin_cos_range( R.xx() ) );
			euler(2) = 0.0;
			euler(3) = 0.0;
		} else if ( R.zz() <= -1 + FLOAT_PRECISION ){
			euler(1) = std::acos( sin_cos_range( R.xx() ) );
			euler(2) = 0.0;
			euler(3) = (T)numeric::constants::d::pi;
		} else {
			T pos_sin_theta = std::sqrt( 1 - R.zz()*R.zz() ); // sin2theta = 1 - cos2theta.
			euler(3) = std::asin( pos_sin_theta );
			if ( R.zz() < 0 ) {
				euler(3) = (T)numeric::constants::d::pi - euler(3);
			}
			euler(1) = std::atan2( R.xz(), -R.yz() );
			euler(2) = std::atan2( R.zx(),  R.zy() );
		}
		euler(1) += euler(1)<0.0 ? (T)numeric::constants::d::pi_2 : 0.0;
		euler(2) += euler(2)<0.0 ? (T)numeric::constants::d::pi_2 : 0.0;
		return euler;
	}

	xyzVector< T > euler_angles_deg() const {
		return (T)numeric::constants::d::radians_to_degrees * euler_angles_rad();
	}
	void from_euler_angles_rad( xyzVector< T > const & euler ) {
		T const phi( euler( 1 ) ), psi( euler( 2 ) ), theta( euler( 3 ) );
		T const cos_phi( std::cos( phi ) ),    sin_phi( std::sin( phi ) );
		T const cos_psi( std::cos( psi ) ),    sin_psi( std::sin( psi ) );
		T const cos_theta( std::cos( theta )), sin_theta( std::sin( theta ) );
		R.xx() =  cos_psi * cos_phi - cos_theta * sin_phi * sin_psi;  R.yx() =  cos_psi * sin_phi + cos_theta * cos_phi * sin_psi;  R.zx() =  sin_psi * sin_theta;
		R.xy() = -sin_psi * cos_phi - cos_theta * sin_phi * cos_psi;  R.yy() = -sin_psi * sin_phi + cos_theta * cos_phi * cos_psi;  R.zy() =  cos_psi * sin_theta;
		R.xz() =                sin_theta * sin_phi;                  R.yz() =                     -sin_theta * cos_phi;            R.zz() =        cos_theta;
	}
	void from_euler_angles_deg( xyzVector< T > const & euler ) {
		xyzVector< T > euler_rad( euler );
		euler_rad *= (T)numeric::constants::d::degrees_to_radians;
		from_euler_angles_rad( euler_rad );
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
	xyzVector<T> const test_point=xyzVector<T>(Real(1.0),Real(3.0),Real(10.0))
){
	xyzTransform<T> I;
	utility::vector1<xyzVector<T> > seenit;
	for(int i0a = 0; i0a < 2; ++i0a){ xyzTransform<T> x0a( (i0a==1 ? G1 : I) *   I );
	for(int i0b = 0; i0b < 2; ++i0b){ xyzTransform<T> x0b( (i0b==1 ? G2 : I) * x0a );
	for(int i1a = 0; i1a < 2; ++i1a){ xyzTransform<T> x1a( (i1a==1 ? G1 : I) * x0b );
	for(int i1b = 0; i1b < 2; ++i1b){ xyzTransform<T> x1b( (i1b==1 ? G2 : I) * x1a );
	for(int i2a = 0; i2a < 2; ++i2a){ xyzTransform<T> x2a( (i2a==1 ? G1 : I) * x1b );
	for(int i2b = 0; i2b < 2; ++i2b){ xyzTransform<T> x2b( (i2b==1 ? G2 : I) * x2a );
	for(int i3a = 0; i3a < 2; ++i3a){ xyzTransform<T> x3a( (i3a==1 ? G1 : I) * x2b );
	for(int i3b = 0; i3b < 2; ++i3b){ xyzTransform<T> x3b( (i3b==1 ? G2 : I) * x3a );
	for(int i4a = 0; i4a < 2; ++i4a){ xyzTransform<T> x4a( (i4a==1 ? G1 : I) * x3b );
	for(int i4b = 0; i4b < 2; ++i4b){ xyzTransform<T> x4b( (i4b==1 ? G2 : I) * x4a );
	for(int i5a = 0; i5a < 2; ++i5a){ xyzTransform<T> x5a( (i5a==1 ? G1 : I) * x4b );
	for(int i5b = 0; i5b < 2; ++i5b){ xyzTransform<T> x5b( (i5b==1 ? G2 : I) * x5a );
	for(int i6a = 0; i6a < 2; ++i6a){ xyzTransform<T> x6a( (i6a==1 ? G1 : I) * x5b );
	for(int i6b = 0; i6b < 2; ++i6b){ xyzTransform<T> x6b( (i6b==1 ? G2 : I) * x6a );
	for(int i7a = 0; i7a < 2; ++i7a){ xyzTransform<T> x7a( (i7a==1 ? G1 : I) * x6b );
	for(int i7b = 0; i7b < 2; ++i7b){ xyzTransform<T> x7b( (i7b==1 ? G2 : I) * x7a );
	for(int i8a = 0; i8a < 2; ++i8a){ xyzTransform<T> x8a( (i8a==1 ? G1 : I) * x7b );
	for(int i8b = 0; i8b < 2; ++i8b){ xyzTransform<T> x8b( (i8b==1 ? G2 : I) * x8a );
	for(int i9a = 0; i9a < 2; ++i9a){ xyzTransform<T> x9a( (i9a==1 ? G1 : I) * x8b );
	for(int i9b = 0; i9b < 2; ++i9b){ xyzTransform<T> x9b( (i9b==1 ? G2 : I) * x9a );
	for(int iAa = 0; iAa < 2; ++iAa){ xyzTransform<T> xAa( (i7a==1 ? G1 : I) * x9b );
	for(int iAb = 0; iAb < 2; ++iAb){ xyzTransform<T> xAb( (i7b==1 ? G2 : I) * xAa );
	for(int iBa = 0; iBa < 2; ++iBa){ xyzTransform<T> xBa( (i8a==1 ? G1 : I) * xAb );
	for(int iBb = 0; iBb < 2; ++iBb){ xyzTransform<T> xBb( (i8b==1 ? G2 : I) * xBa );
	for(int iCa = 0; iCa < 2; ++iCa){ xyzTransform<T> xCa( (i9a==1 ? G1 : I) * xBb );
	for(int iCb = 0; iCb < 2; ++iCb){ xyzTransform<T> xCb( (i9b==1 ? G2 : I) * xCa );
		xyzTransform<T> const & x(xCb);
		if(x.t.length() > r) continue;
		xyzVector<T> test = x*test_point;
		bool redundant(false);
		for(typename utility::vector1<xyzVector<T> >::const_iterator i = seenit.begin(); i != seenit.end(); ++i){
			if( i->distance_squared(test) < 0.001 ) redundant = true;
		}
		if(redundant) continue;
		seenit.push_back(test);
		*container++ = x;
	}} if( 1==N) return;
	}} if( 2==N) return;
	}} if( 3==N) return;
	}} if( 4==N) return;
	}} if( 5==N) return;
	}} if( 6==N) return;
	}} if( 7==N) return;
	}} if( 8==N) return;
	}} if( 9==N) return;
	}} if(10==N) return;
	}} if(11==N) return;
	}} if(12==N) return;
	}}
}

template<typename T, class OutputIterator>
void
expand_xforms(
	OutputIterator container,
	xyzTransform<T> const & G1,
	xyzTransform<T> const & G2,
	int N=5,
	Real r=9e9,
	xyzVector<T> const test_point=xyzVector<T>(Real(1.0),Real(3.0),Real(10.0))
){
	expand_xforms(container,G1,G2,xyzTransform<T>(),N,r,test_point);
}



} // namespace numeric


#endif // INCLUDED_numeric_xyzTransform_HH
