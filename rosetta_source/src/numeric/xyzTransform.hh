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
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>

// C++ headers
#include <cassert>


namespace numeric {

template< typename T >
class xyzTransform {
public:
	// it would be insanity to implement this differently, so these are public members

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
	static inline    X rot       ( MCR rot,          VCR cen){ return X( rot, rot*-cen + cen );}
	static inline    X rot       ( VCR axs, TCR ang, VCR cen){ return rot(rotation_matrix        (axs,ang),cen); }
	static inline    X rot_deg   ( VCR axs, TCR ang, VCR cen){ return rot(rotation_matrix_degrees(axs,ang),cen); }

	inline V     xform( VCR v ) const { return R*v+t; }
	inline V inv_xform( VCR v ) const { return R.transposed()*(v-t); }
};




} // namespace numeric


#endif // INCLUDED_numeric_xyzTransform_HH
