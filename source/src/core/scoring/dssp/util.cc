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
/// @details
/// @author Oliver Lange


// Unit Headers
#include <core/scoring/dssp/util.hh>

// Project Headers
#include <core/pose/Pose.hh>


// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>

//numeric headers

//// C++ headers
// #include <string>

#include <core/chemical/ResidueType.hh>
#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace dssp {

using namespace core;
using namespace pose;
using namespace kinematics;

static basic::Tracer tr( "core.scoring.dssp" );

void
get_CA_vectors(
	PointList const & ca1, // pass by reference, so no tricks:: 3x3
	PointList const & ca2, // pass by reference, so no tricks:: 3x3
	Vector & a,
	Vector & b,
	Vector & c
)
{
	/*       a goes from c-alpha #1 to c-alpha #3 */
	a = ca1[3]-ca1[1];
	a.normalize();

	/*       b gives direction of pleat for ca1 c-alphas */
	Vector b1 = ca1[2] - ca1[1];
	Vector b2 = ca1[2] - ca1[3];
	b = b1 + b2 ;
	b.normalize();

	/*       c goes from ca1 triple to ca2 triple (central alpha-carbons) */
	c = ca2[2] - ca1[2];
	c.normalize();
}


void
get_pairing_geometry(
	pose::Pose const& pose,
	Size const res1,
	Size const res2,
	Real& orientation,
	Real& pleating1,
	Real& pleating2
)
{
	runtime_assert( res1>1 && res1 < pose.size() &&
		res2 > res1 && res2 < pose.size() );

	chemical::ResidueType const& rt1 ( pose.residue_type ( res1 ) );
	chemical::ResidueType const& rt2 ( pose.residue_type ( res2 ) );

	PointList pCA1(3); //space for 3 CAs positions
	PointList pCA2(3);

	// read CAs of 3 consecutive residues
	for ( int i = -1 ; i<=1 ; i++ ) {
		id::AtomID CA1( rt1.atom_index ("CA") , res1+i );
		id::AtomID CA2( rt2.atom_index ("CA") , res2+i );
		pCA1[ i+2 ] = pose.xyz( CA1 );
		pCA2[ i+2 ] = pose.xyz( CA2 );
	};

	Vector dCaCa = pCA1[ 2 ] - pCA2[ 2 ];
	if ( dCaCa.length() > 10.5 ) {
		tr.Warning << "the CA-CA distance for pairing " << res1 << " " << res2 << " is " << dCaCa.length() << std::endl;
		// I put this in because I had this case once, and such an exit would have saved me time.
		//If you have longer distances you probably don't want to choose pleating and orientation based on the native structure
		// --- which is a temporary convenience hack anyway.
		//  utility_exit_with_message(" the CA-CA distance for the chosen pairing is more than 10.5 A check your pairings!!! ");
	}


	Vector a1,b1,c1;
	Vector a2,b2,c2;
	get_CA_vectors( pCA1, pCA2, a1,b1,c1 );
	get_CA_vectors( pCA2, pCA1, a2,b2,c2 );

	orientation = dot_product( a1, a2 );

	Vector ab1,ab2;
	cross(a1,b1,ab1);
	cross(a2,b2,ab2);

	Real const d1 = dot_product(ab1,c1);
	Real const d2 = dot_product(ab2,c2);

	pleating1 = d1;

	if ( orientation < 0 ) {
		pleating2 =  d2; // antiparallel
	} else {
		pleating2 = -d2;
	}
}


void
get_pleating(
	pose::Pose const& pose,
	Size const pos1,
	Size const pos2,
	Size &orientation,
	Size &pleating
)
{

	//Why did this have to get so complicated?
	// Its actually a pretty simple concept!
	//
	// For some reasons, get_pairing_geometry flips
	// pleating2 depending on the orientation --
	// in ideal strand pairs, pleating1 and pleating2 then have the same sign.
	//
	// But for some twisted strand pairs (see, e.g, 22,48 in 1brs.pdb),
	// the numbers get a little crazy...


	Real forientation, pleating1, pleating2;
	if ( pos1 < pos2 ) {
		get_pairing_geometry( pose, pos1, pos2, forientation,pleating1,pleating2);

		//This isn't always quite true...
		//  runtime_assert( pleating1 * pleating2 > 0.0 );
		pleating = ( (pleating1+pleating2) < 0 ? 1 : 2 );
	} else {
		get_pairing_geometry( pose, pos2, pos1, forientation,pleating1,pleating2);


		//This isn't always quite true...
		//  runtime_assert( pleating1 * pleating2 > 0.0 );

		if ( forientation < 0 ) {
			// pleating for anti-parallel pairings is preserved when we
			// interchange positions
			pleating = ( (pleating1+pleating2) < 0 ? 1 : 2 );
		} else {
			// pleating for parallel pairings is reversed when we
			// interchange positions
			pleating = ( (pleating1+pleating2) < 0 ? 2 : 1 );
		}
	}
	tr.Debug << " orientation " << forientation << " pleating " << pleating << std::endl;
	orientation = forientation < 0 ? 1 : 2;
}


}
}
}
