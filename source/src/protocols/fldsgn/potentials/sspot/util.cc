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
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )

#include <protocols/fldsgn/potentials/sspot/util.hh>
#include <core/pose/Pose.hh>

// numeric headers
#include <numeric/conversions.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyzVector.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>


namespace protocols {
namespace fldsgn {
namespace potentials {
namespace sspot {

using Size = core::Size;
using Real = core::Real;
using Vector = core::Vector;
using Pose = core::pose::Pose;

/// @details define a coordinate system with Z axis along cen1->a2 (v1),
/// xz plane defined by cen1->a2 and v21. and origin at cen1.
/// Find the spherical coordinates phi,and theta of point a4 if
/// vector cen2->a4 was moved such that cen2=cen1
/// @param[out] v21  vector connecting midpoints
void
spherical(
	Vector const & a2,
	Vector const & a4,
	Real & phi,
	Real & theta,
	Vector const & cen1,
	Vector const & cen2,
	Vector const & v21
)
{
	using numeric::conversions::degrees;
	using numeric::sin_cos_range;

	Vector v1( a2 - cen1 );   // v1 vector from center to end of dimer1 vector
	Vector v2( a4 - cen2 );   // v2 vector from center to end of dimer2 vector

	Vector const uz( v1.normalized_or_zero() );              // unit vector along v1 = uz
	Vector const uy( uz.cross( v21 ).normalized_or_zero() ); // unit vector perpendicular v21 v1 plane
	Vector const ux( uy.cross( uz ).normalized_or_zero() );  // unit vector to define coordinate system

	Real const v2x = v2.dot( ux ); //v2(1)*ux(1) + v2(2)*ux(2) + v2(3)*ux(3); // v2x=v2.ux
	Real const v2y = v2.dot( uy ); //v2(1)*uy(1) + v2(2)*uy(2) + v2(3)*uy(3); // v2y=v2.uy
	Real const v2z = v2.dot( uz ); //v2(1)*uz(1) + v2(2)*uz(2) + v2(3)*uz(3); // v2z=v2.uz

	Vector const u21( v21.normalized_or_zero() ); // unit vector along v21

	phi = 200.0; // why set to 200?
	//  if v2y = 0, v2 lies in xz plane and phi is zero; if v2x=0, v2 lies in yz plane and phi is 90 or -90
	if ( v2y != 0.0 && v2x != 0.0 ) {
		phi = degrees( std::atan2( v2y, v2x ) );
	}
	if ( phi > 180.0 ) {
		phi -= 360.0f;
	} else if ( phi < -180.0 ) {
		phi += 360.0f;
	}

	Real r1 = v2.length();
	Real const v2z_r1 = v2z/r1;
	// AMW: cppcheck correctly notes that either the first part of this check is useless
	// or we might have just divided by zero
	if ( r1 != 0.0 && std::abs( v2z_r1 ) <= 1.0 ) {
		theta = degrees( numeric::arccos( v2z_r1 ) ); // std::acos( sin_cos_range( v2z_r1 ) ) );
	} else if ( v2z_r1 > 1.0 ) {
		theta = 0.0;
	} else {
		theta = 180.0;
	}

}


/// @brief identifies the sequence separation along the fold tree
/// lin the old way to calculate the sequence separation takes an asumption of no-break chain
/// lin when there is chain break between pos1 and pos2, we add a gap to make a correct calculation in ss energy
Size
get_foldtree_seqsep( Pose const & pose, Size pos1, Size pos2, Size gap_size )
{
	if ( pose.fold_tree().is_simple_tree() ) return std::abs( int( pos1 ) - int( pos2 ) );

	Size begin ( std::min(pos1,pos2) );
	Size end   ( std::max(pos1,pos2) );
	bool is_break ( false );

	for ( Size i = begin; i < end; ++i ) {
		//if( pose.fold_tree().is_cutpoint(i) ) { is_break=true; break; }
		if ( pose.residue_type(i).is_terminus() ) { is_break=true; break; }
	}
	if ( is_break ) {
		return end - begin + gap_size;
	} else {
		return end - begin;
	}
}


} // ns sspot
} // ns potentials
} // ns fldsgn
} // ns protocols

