// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/scoring/sasa.cc
/// @brief  routines which calculate solvent accessible surface area
/// @author Jeff Gray
/// @author Ron Jacak (hydrophobic sasa method, comments, unit tests)

#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <basic/database/open.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
// AUTO-REMOVED #include <core/pose/PDBInfo.hh> // temp
#include <core/scoring/sasa.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <basic/Tracer.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/ubyte.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/Fmath.hh>
#include <ObjexxFCL/format.hh>

// Numeric Headers
#include <numeric/constants.hh>
#include <numeric/trig.functions.hh>

// Utility Headers
#include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/string_util.hh> // temp

// C++ Headers
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>

//#define FILE_DEBUG 1

//#ifdef FILE_DEBUG
static thread_local basic::Tracer TR( "core.scoring.sasa" );
//#endif

using namespace ObjexxFCL::format;

namespace core {
namespace scoring {


// this lookup table is used in sasa computation (also in void.cc)
// the first 8 ints represent the number of 1 bits you would find in the binary equivalents of the decimal numbers 0-7.
// e.g. 1 in binary is 0000 0001, so 1 '1'-bit
//      7 in binary is 0000 0111, so 3 '1'-bits
short const bit_count[] = { // lookup table for number of 1 bits in a ubyte
	0,1,1,2,1,2,2,3, 1,2,2,3,2,3,3,4,   1,2,2,3,2,3,3,4, 2,3,3,4,3,4,4,5,  // 0x 1x
	1,2,2,3,2,3,3,4, 2,3,3,4,3,4,4,5,   2,3,3,4,3,4,4,5, 3,4,4,5,4,5,5,6,  // 2x 3x
	1,2,2,3,2,3,3,4, 2,3,3,4,3,4,4,5,   2,3,3,4,3,4,4,5, 3,4,4,5,4,5,5,6,  // 4x 5x
	2,3,3,4,3,4,4,5, 3,4,4,5,4,5,5,6,   3,4,4,5,4,5,5,6, 4,5,5,6,5,6,6,7,  // 6x 7x
	1,2,2,3,2,3,3,4, 2,3,3,4,3,4,4,5,   2,3,3,4,3,4,4,5, 3,4,4,5,4,5,5,6,  // 8x 9x
	2,3,3,4,3,4,4,5, 3,4,4,5,4,5,5,6,   3,4,4,5,4,5,5,6, 4,5,5,6,5,6,6,7,  // Ax Bx
	2,3,3,4,3,4,4,5, 3,4,4,5,4,5,5,6,   3,4,4,5,4,5,5,6, 4,5,5,6,5,6,6,7,  // Cx Dx
	3,4,4,5,4,5,5,6, 4,5,5,6,5,6,6,7,   4,5,5,6,5,6,6,7, 5,6,6,7,6,7,7,8,  // Ex Fx
};


int const num_bytes = 21;
int const num_phi = 64;
int const num_theta = 64;
int const num_overlaps = 100;
int const num_orientations = 162;
int const maskbits = 162;

// This database file is used to find the number of the closest point to the point of intersection between two spheres.
// The surface of each atom is covered by 162 points. Very rarely is it the case that two spheres will intersect exactly
// on a point.  This database table is used to map the actual point of intersection to the closest predetermined dot.
ObjexxFCL::FArray2D_int angles( num_phi, num_theta );

// a 2D FArray of ubytes, that is num_bytes (21) by overlaps * orientations.
// Each line of this database file contains 21 ints ranging from 0-255. These values in binary represent the mask
// that should get applied to the state of an atoms current set of dots. So if the line is 21 0's that means that no
// bits in the atoms state should be changed. (No overlap occurs between these two atoms.)  On the other hand, if the
// line is 21 255's, that means the atom is completely covered by the other atom.
ObjexxFCL::FArray2D_ubyte masks( num_bytes, num_overlaps * num_orientations );


///
/// @begin sasa.cc::input_sasa_dats
///
/// @brief
/// Reads in the SASA database files sampling/SASA-angles.dat and sampling/SASA-masks.dat into FArrays above.
///
void input_sasa_dats() {

	// this booleans checks whether we've done this already or not. if we have, then return immediately.
	static bool init = false;
	if ( init )
		return;
	init = true;

	//j inputting the masks. they are 21 ubytes long, 162x100 (see header). expects file to be complete
	utility::io::izstream masks_stream( basic::database::full_name("sampling/SASA-masks.dat" ) );

	//ronj iterate over all 16200 bit masks in the database file and save them to the FArray 'masks'
	//ronj the lines have 21 space-delimited DECIMAL values which we need to convert to binary values that are 8-bits long
	//ronj so save the int (or short) somewhere and then cast that value into an unsigned char

	for ( int ii=1; ii <= num_overlaps * num_orientations; ++ii ) {
		short tmp; // a 'short' is a int type that is usually 16-bits
		for ( int jj = 1; jj <= num_bytes; ++jj ) {
			//ronj operator>> of izstream (aka the extraction operator) performs an input operation on a stream generally
			//ronj involving some interpretation of the data
			masks_stream >> tmp;
			masks(jj,ii) = static_cast< unsigned char >( tmp );
		}
		masks_stream >> skip; // calls skip() on the izstream which skips the rest of the line and the line terminator
	}
	masks_stream.close();

	//ronj iterate over all 64 lines of the sampling/SASA-angles.database file and save them to the FArray 'angles'
	//ronj these lines have 64 space-delimited ints which need to be stored
	utility::io::izstream angles_stream( basic::database::full_name( "sampling/SASA-angles.dat" ) );

	for ( int ii = 1; ii <= num_phi; ++ii ) {
		for ( int jj = 1; jj <= num_theta; ++jj ) {
			angles_stream >> angles( ii, jj );
		}
		angles_stream >> skip;
	}
	angles_stream.close();




}


///
/// @begin sasa.cc::get_overlap
///
/// @detailed
/// getting overlap from a to b (or i to j, as the atoms are referred to in calc_per_atom_sasa below).
/// this returns the degree of overlap between two atoms adapted from erics code in area.c GetD2 and returns value
/// from 1 to 100. This calculation is based on the law of cosines.
/// See LeGrand and Merz, Journal of Computational Chemistry 14(3):349-52 (1993).
/// Note that equation (4) is wrong, the denominator should be 2*ri*riq  instead of  2*ri*rq   (j)
///
/// The function gets passed in the sasa radius of atom i (plus the probe radius), the sasa radius of atom j (plus
/// the probe radius), the distance between the atom centers, and a reference to the degree of overlap (represented
/// as an int). The degree of overlap that's returned can be thought of as how much of atom a is covered by atom b.
/// A value of 100 means that atom a is completely covered up by atom b. A value of 1 means that not much of the surface
/// of 'a' is covered up by 'b'.
/// The law of cosines relates the cosine of one angle of a triangle to the lengths of its sides. More specifically,
/// c^2 = a^2 + b^2 - 2*a*b*cos theta, where theta is the angle between sides a and b.  For the function we want to
/// compute the angle of the cone of intersection between spheres 'a' and 'b'.  Let the radius of atom a be ri, and the
/// radius of atom b be rq, and the distance between atom centers be riq.  Let the angle between ri and riq be theta_iq.
/// The cosine of theta_iq will be equivalent to ( ri^2 + riq^2 - rq^2 ) / 2 * ri * riq
///
void
get_overlap( Real const radius_a, Real const radius_b, Real const distance_ijxyz, int & degree_of_overlap ) {

	//j min distance cutoff
	Real epsilon = 0.01;

	if ( distance_ijxyz < epsilon ) {
		//j atoms too close, causes round off error. use this cutoff.
		//ronj not sure what the rationale for using the cutoffs below is.
		if ( radius_a < radius_b ) {
			degree_of_overlap = 100;
		} else {
			degree_of_overlap = 1;
		}

	} else if ( radius_b + distance_ijxyz <= radius_a ) {
		//j If atom a completely engulfs atom b, consider a to have no overlap due to atom b.
		degree_of_overlap = 1;

	} else if ( radius_a + distance_ijxyz <= radius_b ) {
		//j If atom a is completely engulfed by atom b, then turn it completely off (i.e. d2 = 99).
		degree_of_overlap = 100;

	} else {
		//j Otherwise, compute the amount of overlap using the law of cosines. "costh" is the angle of the cone of
		//j intersection that atom b imposes on atom a.
		Real cosine_theta_iq = ( ( radius_a*radius_a ) + ( distance_ijxyz*distance_ijxyz ) - ( radius_b*radius_b ) ) / ( 2 * radius_a * distance_ijxyz );

		//ronj when 'b' and 'a' are barely overlapping, the theta angle will be small. as 'b' moves closer to 'a', the point of
		//ronj intersection moves farther and farther away from 'b' causing the theta angle to become wider and wider.
		//ronj the theta angle will vary from 0deg to 180deg. the cos of theta will vary from 1 to 0 to -1.
		//ronj therefore, subtracting cos theta from 1 and multiplying by 50 will result in the degree_of_overlap varying
		//ronj from 1 (when cos theta is 1, or angle is 0) to 100 (when cos theta is -1, or theta is 180deg).
		degree_of_overlap = static_cast< int >( (1.0f - cosine_theta_iq) * 50 ) + 1;

		if ( degree_of_overlap > 100 ) {
			degree_of_overlap = 100;

		} else if ( degree_of_overlap < 0 ) {
			//j We already hopefully accounted for this possibility by requiring that distance_ijxyz > epsilon,
			//j but in case not we don't want a potential bug to go unnoticed.
			std::cout << "Problem in calculating overlap between atoms. Terminating calculation.\n";
			std::cout << "radius_a: " << SS( radius_a ) << ", radius_b: " << SS( radius_b )
				<< ", distance_ijxyz: " << SS( distance_ijxyz ) << ", cosine theta: " << SS( cosine_theta_iq ) << std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}
	}
}


///
/// @begin sasa.cc::get_orientation
///
/// @brief
/// Gets the orientation of a to b (i to j, see below). Does this by calculating two angles, aphi and theta. (j)
///
/// @detailed
/// ronj This function is used to get two indexes (phi and theta) which are used to get the index of a dot on the
/// ronj surface of the 'a' sphere. When calculating how much surface area sphere b covers on a, we can get the degree
/// ronj of overlap from the function above, but it's not necessarily the case that the vector that connects the center
/// ronj of atom 'a' and atom 'b' goes through one of the predetermined dot locations on the surface of 'a'. In fact,
/// ronj it's very unlikely that the vector goes through a predetermined dot.  Instead, what is done is the actual point
/// ronj of intersection (the outermost point of a on the line from the center of 'a' to center of 'b') is converted
/// ronj to spherical polar coordinates. Then, the values are used to find the location of the closest predetermined
/// ronj point on the surface of 'a' using a lookup table. So what this function needs to do is convert the
/// ronj cartesian coordinate of the actual point of intersection into polar coordinates.
/// ronj
/// ronj To get the spherical, polar coordinates of a cartesian point x,y,z, use these equations:
/// ronj r = sqrt( x^2 + y^2 + z^2 )
/// ronj theta = arccos( z / r )
/// ronj phi = arctan( y / x )
///
/// ronj Then, once we have the true phi and theta, we need to translate this into an index (or offset) for the correct
/// ronj value in the database file. There are 64 phi angle bin and 64 theta bins in the database file sampling/SASA-angles.dat.
/// ronj We need to convert the phi and theta into indexes for this file by multiplying them by num_phi / 2*pi.
//  ronj pi_2 is 2*pi
///
/// ronj Note: I think phi and theta have been reversed in the function below. The code below uses the following:
/// ronj phi = arccos( z )
/// ronj theta = arctan( y / x )
///
/// ronj After a couple of weeks trying to write tests for this function, I have been unsuccessful in figuring out why
/// ronj it's doing what it does.  Despite using the wrong equations, it seems to work.  Comparing the total residue
/// ronj SASA values calculated by calc_per_atom_sasa() below results in a correlation of 0.98 against what the program
/// ronj NACCESS finds for the same residues. This test was done on a small 110aa protein.  I also looked at the per-atom
/// ronj total SASA and the correlation for all atoms (mini v. NACCESS) was approximately 0.94. I'm using exactly the same
/// ronj van der Waals radii for both programs so I feel like the correlations should be 1.0. Explanations for the
/// ronj differences can be 1) this method is doing something wrong in calculating the closest surface point, 2) this
/// ronj method is correct but the masks that are in the database are not aligned to the surface points correctly, 3) the
/// ronj differences are solely due to the different way that the two program calculate surface area.
///
///
void get_orientation( Vector const & a_xyz, Vector const & b_xyz, int & phi_index, int & theta_index, Real distance_ijxyz ) {

	using namespace numeric::constants::d;
	using numeric::sin_cos_range;

	//static FArray1D_float diff( 3 ); //apl allocate this once only
	//pb -- static can cause problems in multi-threading
	Vector ab_diff( ( a_xyz - b_xyz ) / distance_ijxyz );

	//j now figure out polar values of phi and theta. Normalize the difference. first get the length of the vector.
	//vector_normalize(ab_diff);

	//j figuring phi
	//ronj sin_cos_range adjust the range of the passed in argument to a valid [-1,1]. acos(x) returns the arccos of x
	//ronj expressed in radians. (arccos is the inverse operation of cos.) the return value is in the range [0,pi].
	Real phi = std::acos( sin_cos_range( ab_diff(3) ) );
#ifdef FILE_DEBUG
	Real temp_phi = phi;
#endif
	//ronj if z = -1, then phi becomes arccos( -1 ) = pi
	//ronj if z = -0.5, then phi becomes arccos( -0.5 ) = 2pi/3
	//ronj if z = 0, then phi becomes arccos( 0 ) = pi/2
	//ronj if z = 0.5, then phi becomes arccos( 0.5 ) = pi/3
	//ronj if z = 1, then phi becomes arccos( 1 ) = 0

	//ronj phi ranges from [0,pi], so...
	//ronj if phi = 0, then p becomes 0 * ( num_phi / 2*pi ) = 0
	//ronj if phi = pi/2, then p becomes pi/2 * ( num_phi / 2*pi ) = 64 / 4 = 16
	//ronj if phi = pi, then p becomes pi * ( num_phi / 2*pi ) = num_phi / 2 = 32
	phi = phi * ( num_phi / pi_2 );
#ifdef FILE_DEBUG
	Real precasted_phi = phi;
#endif
	phi_index = static_cast< int >( phi );
	++phi_index; // for fortran goes from 1 to n
	if ( phi_index > num_phi )
		phi_index = 1; // reset back to index of 1?

	//j figuring theta
	//ronj atan2 is the arc tangent of y/x, expressed in radians.
	//ronj atan2 ranges from [-pi,pi], so...
	//ronj if atan2 = -pi, then t becomes -pi * ( num_theta / 2*pi ) = -64 / 2 = -32 += 64 = 32
	//ronj if atan2 = -pi/2, then t becomes -pi/2 * ( num_theta / 2*pi ) = -64 / 4 = -16 += 64 = 48
	//ronj if atan2 = 0, then t becomes 0 * ( num_theta / 2*pi ) = 0
	//ronj if atan2 = pi/2, then t becomes pi/2 * ( num_theta / 2*pi ) = 64 / 4 = 16
	//ronj if atan2 = pi, then t becomes pi * ( num_theta / 2*pi ) = 64 / 2 = 32
	Real theta = std::atan2( ab_diff(2), ab_diff(1) );
#ifdef FILE_DEBUG
	Real temp_theta = theta;
#endif

	theta = theta * ( num_theta / pi_2 );
#ifdef FILE_DEBUG
	Real precasted_theta = theta;
#endif
	theta_index = static_cast< int >( theta );
	++theta_index; // for fortran goes from 1 to n
	if ( theta_index < 1 ) { // as in the case when atan2 returns a negative...
		theta_index += num_theta; // let t = -32; -32 + 64 = 32
	} else if ( theta_index > num_theta ) {
		theta_index = 1; // reset back to index of 1?
	}

#ifdef FILE_DEBUG
	/*TR << "get_orientation(): ";
	TR << "a: " << F(5,2,a_xyz.x()) << "," << F(5,2,a_xyz.y()) << "," << F(5,2,a_xyz.z())
		<< ", b: " << F(5,2,b_xyz.x()) << "," << F(5,2,b_xyz.y()) << "," << F(5,2,b_xyz.z())
		<< ", a-b: " << F(6,3,ab_diff.x()) << "," << F(6,3,ab_diff.y()) << "," << F(6,3,ab_diff.z())
		<< ", distance_ijxyz: " << F( 5,3,distance_ijxyz )
		<< ", phi: " << temp_phi << ", unrounded_phi: " << precasted_phi << ", phi_index: " << SS( phi_index )
		 << ", theta: " << temp_theta << ", unrounded_theta: " << precasted_theta << ", theta_index: " << SS( theta_index )
		 << std::endl;
	TR << "get_orientation(): closest dot: " << angles( phi_index, theta_index ) << std::endl;*/
#endif
}


///
/// @begin sasa.cc::get_2way_orientation
///
/// @brief
/// Gets the orientation of a to b (i to j, see below). Does this by calculating two angles, aphi and theta. (j)
///
/// @detailed
/// ronj This function is the same as the function above but get the orientation of a to b simultaneously with the
/// ronj orientation of b to a.  The same result could be achieved by making two separate get_2way_orientation() calls
/// ronj but this method does it more efficiently by avoiding an atan2 and acos call.  Instead, once you compute the
/// ronj phi and theta for a on b, you can add/subtrate pi factors to get the phi and theta for b on a.
/// ronj Still not sure how this method returns the correct values, though.
///
void
get_2way_orientation( Vector const & a_xyz, Vector const & b_xyz, int & phi_a2b_index, int & theta_a2b_index, int & phi_b2a_index, int & theta_b2a_index, Real distance_ijxyz ) {

	using namespace numeric::constants::d;
	using numeric::sin_cos_range;

	//static FArray1D_float diff( 3 ); //apl allocate this once only
	//pb -- static can cause problems in multi-threading
	Vector ab_diff( ( a_xyz - b_xyz ) / distance_ijxyz );

	//j now figure out polar values of phi and theta. Normalize the difference. first get the length of the vector.
	//vector_normalize(ab_diff);

	//j figuring phi, for a2b
	//ronj sin_cos_range adjust the range of the passed in argument to a valid [-1,1]. acos(x) returns the arccos of x
	//ronj expressed in radians. (arccos is the inverse operation of cos.) the return value is in the range [0,pi].
	Real phi_a2b = std::acos( sin_cos_range( ab_diff(3) ) );
#ifdef FILE_DEBUG
	Real temp_phi_a2b = phi_a2b;
#endif
	//ronj if z = -1, then phi_a2b becomes arccos( -1 ) = pi
	//ronj if z = -0.5, then phi_a2b becomes arccos( -0.5 ) = 2pi/3
	//ronj if z = 0, then phi_a2b becomes arccos( 0 ) = pi/2
	//ronj if z = 0.5, then phi_a2b becomes arccos( 0.5 ) = pi/3
	//ronj if z = 1, then phi_a2b becomes arccos( 1 ) = 0

	//ronj phi ranges from [0,pi], so...
	//ronj if phi_a2b = 0, then p becomes 0 * ( num_phi / 2*pi ) = 0
	//ronj if phi_a2b = pi/2, then p becomes pi/2 * ( num_phi / 2*pi ) = 64 / 4 = 16
	//ronj if phi_a2b = pi, then p becomes pi * ( num_phi / 2*pi ) = num_phi / 2 = 32
	phi_a2b = phi_a2b * ( num_phi / pi_2 );
#ifdef FILE_DEBUG
	Real precasted_phi_a2b = phi_a2b;
#endif
	phi_a2b_index = static_cast< int >( phi_a2b );
	++phi_a2b_index; // for fortran goes from 1 to n
	if ( phi_a2b_index > num_phi )
		phi_a2b_index = 1; // reset back to index of 1?

	// can take a shortcut to get phi_b2a
	Real phi_b2a = num_phi / 2 - phi_a2b;
#ifdef FILE_DEBUG
	Real temp_phi_b2a = phi_b2a;
#endif
	phi_b2a_index = static_cast< int >( phi_b2a );
	++phi_b2a_index; // for fortran goes from 1 to n
	if ( phi_b2a_index > num_phi )
		phi_b2a_index = 1;


	//j figuring theta
	//ronj atan2 is the arc tangent of y/x, expressed in radians.
	//ronj atan2 ranges from [-pi,pi], so...
	//ronj if atan2 = -pi, then t becomes -pi * ( num_theta / 2*pi ) = -64 / 2 = -32 += 64 = 32
	//ronj if atan2 = -pi/2, then t becomes -pi/2 * ( num_theta / 2*pi ) = -64 / 4 = -16 += 64 = 48
	//ronj if atan2 = 0, then t becomes 0 * ( num_theta / 2*pi ) = 0
	//ronj if atan2 = pi/2, then t becomes pi/2 * ( num_theta / 2*pi ) = 64 / 4 = 16
	//ronj if atan2 = pi, then t becomes pi * ( num_theta / 2*pi ) = 64 / 2 = 32
	Real theta_a2b = std::atan2( ab_diff(2), ab_diff(1) );
#ifdef FILE_DEBUG
	Real temp_theta_a2b = theta_a2b;
#endif

	theta_a2b = theta_a2b * ( num_theta / pi_2 );
#ifdef FILE_DEBUG
	Real precasted_theta_a2b = theta_a2b;
#endif
	theta_a2b_index = static_cast< int >( theta_a2b );
	++theta_a2b_index; // for fortran goes from 1 to n
	if ( theta_a2b_index < 1 ) { // as in the case when atan2 returns a negative...
		theta_a2b_index += num_theta; // let t = -32; -32 + 64 = 32
	} else if ( theta_a2b_index > num_theta ) {
		theta_a2b_index = 1; // reset back to index of 1?
	}

	// can use a shortcut to get theta_b2a
	Real theta_b2a = num_theta / 2.0f + theta_a2b;
	if ( theta_b2a > num_theta / 2.0f )
		theta_b2a -= num_theta;
#ifdef FILE_DEBUG
	Real temp_theta_b2a = theta_b2a;
#endif

	theta_b2a_index = static_cast< int >( theta_b2a );
	++theta_b2a_index; // for fortran goes from 1 to n
	if ( theta_b2a_index < 1 ) {
		theta_b2a_index += num_theta;
	} else if ( theta_b2a_index > num_theta ) {
		theta_b2a_index = 1;
	}


#ifdef FILE_DEBUG
	/*TR << "get_2way_orientation(): ";
	TR << "a: " << F(5,2,a_xyz.x()) << "," << F(5,2,a_xyz.y()) << "," << F(5,2,a_xyz.z())
		<< ", b: " << F(5,2,b_xyz.x()) << "," << F(5,2,b_xyz.y()) << "," << F(5,2,b_xyz.z())
		<< ", a-b: " << F(6,3,ab_diff.x()) << "," << F(6,3,ab_diff.y()) << "," << F(6,3,ab_diff.z())
		<< ", distance_ijxyz: " << F( 5,3,distance_ijxyz )
		<< ", phi_a2b: " << temp_phi_a2b << ", unrounded_phi_a2b: " << precasted_phi_a2b << ", phi_a2b_index: " << SS( phi_a2b_index )
		 << ", theta_a2b: " << temp_theta_a2b << ", unrounded_theta_a2b: " << precasted_theta_a2b << ", theta_a2b_index: " << SS( theta_a2b_index )
		 << ", phi_b2a: " << temp_phi_b2a << ", phi_b2a_index: " << SS( phi_b2a_index )
		 << ", theta_b2a: " << temp_theta_b2a << ", theta_b2a_index: " << SS( theta_b2a_index )
		 << std::endl;*/
#endif

}


Real
calc_total_sasa( pose::Pose const & pose, Real const probe_radius ) {

	id::AtomID_Map< Real > atom_sasa;
	utility::vector1< Real > rsd_sasa;

	return calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius );
}


Real
calc_per_atom_sasa( pose::Pose const & pose, id::AtomID_Map< Real > & atom_sasa, utility::vector1< Real > & rsd_sasa,
	Real const probe_radius, bool const use_big_polar_H /* = false */ ) {

	id::AtomID_Map< bool > atom_subset;
	atom_subset.clear();
	core::pose::initialize_atomid_map( atom_subset, pose, true ); // jk use all atoms if atom_subset is not specified

	return calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius, use_big_polar_H, atom_subset );
}

Real
calc_per_atom_sasa_sc( pose::Pose const & pose, utility::vector1< Real > & rsd_sasa,bool normalize) {
	id::AtomID_Map< bool > atom_subset;
	id::AtomID_Map< Real > atom_sasa;
	static Real const  probe_radius=1.4;
	static bool const use_big_polar_H=false;
	atom_subset.clear();
	core::pose::initialize_atomid_map( atom_subset, pose, true);
	Real total_sasa=calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius, use_big_polar_H, atom_subset,true /*use_naccess_sasa_radii*/ );

	for(Size i=1;i<=atom_sasa.n_residue();++i) {
		rsd_sasa[i]=0;
		conformation::Residue const & rsd( pose.residue(i) );
		for(Size ii=1;ii<=atom_sasa.n_atom(i);++ii) {
			if(!rsd.atom_is_backbone(ii)) { // || (rsd.name1()=='G' && rsd.atom_name(ii).compare(1,2,"CA")==0) ){
				id::AtomID const id(ii,i);
				rsd_sasa[i]+=atom_sasa[ id ];
			} else {
				if(rsd.name1()=='G' && rsd.atom_name(ii).compare(1,2,"CA")==0) {
						id::AtomID const id(ii,i);
						rsd_sasa[i]+=atom_sasa[ id ]*1.67234+0.235839; //Somehow Gly area is too small with *1.67234+0.235839;
					}
			//			std::cout << "ATOM_SASA: " << rsd.name3() << " " << rsd.atom_name(j) << " " << id << " " << atom_sasa[ id ] << std::endl;
			}
		}
		if(normalize) {
			rsd_sasa[i]/=normalizing_area(rsd.name1());
			rsd_sasa[i]*=100;
		}
	}


	return total_sasa;
}
//ATOM S   2  ALA  107.95   0.0  69.41   0.0   0.00   0.0  69.41   0.0  38.54   0.0  71.38   0.0  36.58   0.0
//ATOM S   2  CYS  134.28   0.0  96.75   0.0   0.00   0.0  96.75   0.0  37.53   0.0  97.93   0.0  36.35   0.0
//ATOM S   2  ASP  140.39   0.0  48.00   0.0  54.69   0.0 102.69   0.0  37.70   0.0  49.24   0.0  91.15   0.0
//ATOM S   2  GLU  172.25   0.0  59.10   0.0  75.64   0.0 134.74   0.0  37.51   0.0  60.29   0.0 111.96   0.0
//ATOM S   2  PHE  199.48   0.0 164.11   0.0   0.00   0.0 164.11   0.0  35.37   0.0 165.25   0.0  34.23   0.0
//ATOM S   2  GLY   80.10   0.0  32.33   0.0   0.00   0.0  32.33   0.0  47.77   0.0  37.55   0.0  42.55   0.0
//ATOM S   2  HIS  182.88   0.0  96.01   0.0  51.07   0.0 147.08   0.0  35.80   0.0  97.15   0.0  85.73   0.0
//ATOM S   2  ILE  175.12   0.0 137.96   0.0   0.00   0.0 137.96   0.0  37.16   0.0 139.14   0.0  35.98   0.0
//ATOM S   2  LYS  200.81   0.0 115.38   0.0  47.92   0.0 163.30   0.0  37.51   0.0 116.57   0.0  84.24   0.0
//ATOM S   2  LEU  178.63   0.0 141.12   0.0   0.00   0.0 141.12   0.0  37.51   0.0 142.31   0.0  36.32   0.0
//ATOM S   2  MET  194.15   0.0 156.64   0.0   0.00   0.0 156.64   0.0  37.51   0.0 157.84   0.0  36.32   0.0
//ATOM S   2  ASN  143.94   0.0  44.98   0.0  61.26   0.0 106.24   0.0  37.70   0.0  46.23   0.0  97.72   0.0
//ATOM S   2  PRO  136.13   0.0 119.90   0.0   0.00   0.0 119.90   0.0  16.23   0.0 120.95   0.0  15.19   0.0
//ATOM S   2  GLN  178.50   0.0  51.03   0.0  89.96   0.0 140.99   0.0  37.51   0.0  52.22   0.0 126.28   0.0
//ATOM S   2  ARG  238.76   0.0  76.60   0.0 124.65   0.0 201.25   0.0  37.51   0.0  77.80   0.0 160.97   0.0
//ATOM S   2  SER  116.50   0.0  46.89   0.0  31.22   0.0  78.11   0.0  38.40   0.0  48.55   0.0  67.95   0.0
//ATOM S   2  THR  139.27   0.0  74.54   0.0  27.17   0.0 101.70   0.0  37.57   0.0  75.72   0.0  63.55   0.0
//ATOM S   2  VAL  151.44   0.0 114.28   0.0   0.00   0.0 114.28   0.0  37.16   0.0 115.47   0.0  35.97   0.0
//ATOM S   2  TRP  249.36   0.0 187.67   0.0  23.60   0.0 211.26   0.0  38.10   0.0 189.67   0.0  59.69   0.0
//ATOM S   2  TYR  212.76   0.0 135.35   0.0  42.03   0.0 177.38   0.0  35.38   0.0 136.50   0.0  76.26   0.0

static utility::vector1<Real> init_normalizing_area_sc() {
	static utility::vector1<Real> area_sc(255,0);
	area_sc[65]= 69.41;   // 69.41;     //A
	area_sc[67]= 96.75;   // 96.75;		 //C
	area_sc[68]=102.69;   // 48.00;		 //D
	area_sc[69]=134.74;   // 59.10;		 //E
	area_sc[70]=164.11;   //164.11;		 //F
	area_sc[71]= 32.33;   // 32.33;		 //G
	area_sc[72]=147.08;   // 96.01;		 //H
	area_sc[73]=137.96;   //137.96;		 //I
	area_sc[75]=163.30;   //115.38;		 //K
	area_sc[76]=141.12;   //141.12;		 //L
	area_sc[77]=156.64;   //156.64;		 //M
	area_sc[78]=106.24;   // 44.98;		 //N
	area_sc[80]=119.90;   //119.90;		 //P
	area_sc[81]=140.99;   // 51.03;		 //Q
	area_sc[82]=201.25;   // 76.60;		 //R
	area_sc[83]= 78.11;   // 46.89;		 //S
	area_sc[84]=101.70;   // 74.54;		 //T
	area_sc[86]=114.28;   //114.28;		 //V
	area_sc[87]=211.26;   //187.67;		 //W
	area_sc[89]=177.38;   //135.35;		 //Y
	//	std::cout << "INIT!!!!!!HELLO" << std::endl;
	return area_sc;
}

static utility::vector1<Real> init_normalizing_area_total() {
	static utility::vector1<Real> area_total(255,0);
	area_total[65]=107.95;     //A
	area_total[67]=134.28;		 //C
	area_total[68]=140.39;		 //D
	area_total[69]=172.25;		 //E
	area_total[70]=199.48;		 //F
	area_total[71]= 80.10;		 //G
	area_total[72]=182.88;		 //H
	area_total[73]=175.12;		 //I
	area_total[75]=200.81;		 //K
	area_total[76]=178.63;		 //L
	area_total[77]=194.15;		 //M
	area_total[78]=143.94;		 //N
	area_total[80]=136.13;		 //P
	area_total[81]=178.50;		 //Q
	area_total[82]=238.76;		 //R
	area_total[83]=116.50;		 //S
	area_total[84]=139.27;		 //T
	area_total[86]=151.44;		 //V
	area_total[87]=249.36;		 //W
	area_total[89]=212.76;		 //Y
	return area_total;
}
Real normalizing_area(char const res) {

	int index=static_cast<unsigned int>(static_cast<unsigned char>(res));
	static utility::vector1<Real> const area_sc=init_normalizing_area_sc();
	static utility::vector1<Real> const area_total=init_normalizing_area_total();


	//std::cout << "INDEX: " << res << " " << index << " " << area_sc[index] << " " << area_total[index] << std::endl;;
	return area_sc[index];


}

Real
calc_per_atom_sasa(
	pose::Pose const & pose,
	id::AtomID_Map< Real > & atom_sasa,
	utility::vector1< Real > & rsd_sasa,
	Real const probe_radius,
	bool const use_big_polar_H,
	id::AtomID_Map< bool > & atom_subset,
	bool const use_naccess_sasa_radii /* =false */,
	bool const expand_polar_radii /* =false */,
	Real const polar_expansion_radius /* =1.0 */,
	bool const include_probe_radius_in_atom_radii, /* = true; used in calc of final sasas */
	bool const use_lj_radii /* = false */
) {

	using core::conformation::Residue;
	using core::conformation::Atom;
	using core::id::AtomID;

	if ( pose.total_residue() < 1 )
		return 0.0; // nothing to do

	// read sasa datafiles
	input_sasa_dats();

	Real const big_polar_H_radius( 1.08 ); // increase radius of polar hydrogens, eg when doing unsatisfied donorH check

	//j setup the radii array, indexed by the atom type int. atom index for looking up an extra data type stored in the AtomTypes
	//ronj reads the values out of the database file sasa_radii.txt in the extras folder of atom_type_sets and stores the values
	//ronj for each atom type into the radii array. each index of the radii array corresponds to some atom type.
	core::chemical::AtomTypeSet const & atom_type_set = pose.residue(1).atom_type_set();

	core::Size SASA_RADIUS_INDEX;
	if ( use_naccess_sasa_radii ) {
		SASA_RADIUS_INDEX = atom_type_set.extra_parameter_index( "NACCESS_SASA_RADIUS" );
	} else {
		SASA_RADIUS_INDEX = atom_type_set.extra_parameter_index( "REDUCE_SASA_RADIUS" );
	}

	utility::vector1< Real > radii( atom_type_set.n_atomtypes() );

	for ( core::Size ii=1; ii <= radii.size(); ++ii ) {
		core::chemical::AtomType const & at( atom_type_set[ii] );

		if ( use_lj_radii ) {
			radii[ii] = atom_type_set[ii].lj_radius();
// 			std::cout << "Using LJ radius  "<< radii[ii] << " instead of " << atom_type_set[ii].extra_parameter( SASA_RADIUS_INDEX ) <<
// 				" for " << atom_type_set[ii].name() << std::endl;
			continue; // note continue, no other checks applied
		}

		radii[ii] = atom_type_set[ii].extra_parameter( SASA_RADIUS_INDEX );

		if ( use_big_polar_H && at.is_polar_hydrogen() && big_polar_H_radius > radii[ii] ) {
			std::cout << "Using " << big_polar_H_radius << " instead of " << radii[ii] << " for atomtype " << at.name() << " in sasa calculation!\n";
			radii[ii] = big_polar_H_radius;
		}

		if ( expand_polar_radii && ( at.element() == "N" || at.element() == "O" ) ) {
#ifdef FILE_DEBUG
			std::cout << "Using " << radii[ii] + polar_expansion_radius << " instead of " << radii[ii] << " for atomtype " << at.name() << " in sasa calculation!\n";
#endif
			radii[ii] = radii[ii] + polar_expansion_radius;
		}

	}

	// create an AtomID_Map which will map atom ids to the mask for that atom. right now, just make all atom ids point
	// to a mask of all zeros. later we'll have to set different masks.
	core::id::AtomID_Map< utility::vector1< ObjexxFCL::ubyte > > atom_masks;
	utility::vector1< ObjexxFCL::ubyte > zero_mask( num_bytes, ObjexxFCL::ubyte( 0 ) );
	core::pose::initialize_atomid_map( atom_masks, pose, zero_mask ); // calls AtomID_Map.Pose.hh::initialize()


	// identify the maxium radius we have for all atoms in the pose.  this will be used to set a cutoff distance for use
	// in skipping residue pairs if their nbr atoms are too far apart.
	core::Real max_radius = 0.0;
	for ( Size ii=1; ii <= pose.total_residue(); ++ii ) {
		Residue const & rsd( pose.residue( ii ) );
		// iterate over all the atoms in residue ii and
		for ( Size jj=1; jj <= rsd.natoms(); ++jj ) {
			max_radius = std::max( max_radius, radii[ rsd.atom(jj).type() ] );
			//ronj assert that the sum of the first 8 bits and the last 8 bits (or first byte and last byte) in the vector
			//ronj of masks for this atom id is less than 0.001.  what is the point of this assertion?
			runtime_assert( ( atom_masks[ AtomID(jj,ii) ][1] + atom_masks[ AtomID(jj,ii) ][num_bytes] ) < 1e-3 );
		}
	}

	core::Real cutoff_distance = 0.0;
	cutoff_distance = 2 * ( max_radius + probe_radius );


	//j now do calculations: get the atom_masks by looping over all_atoms x all_atoms
	for ( Size ii=1; ii <= pose.total_residue(); ++ii ) {
		Residue const & irsd( pose.residue( ii ) );

		//ronj for the other 'j' residue, only iterate over residues which have indexes > residue 'i'
		for ( Size jj=ii; jj <= pose.total_residue(); ++jj ) {
			Residue const & jrsd( pose.residue( jj ) );


			calc_atom_masks(
				irsd, jrsd,
				probe_radius, cutoff_distance, radii,
				atom_subset,
				atom_masks );

		} // iia
	} // ir

	//j calculate the residue and atom sasa
	Real total_sasa( 0.0 );

	rsd_sasa.clear();
	rsd_sasa.resize( pose.total_residue(), 0.0 );

	atom_sasa.clear();
	core::pose::initialize_atomid_map( atom_sasa, pose, (Real) -1.0 ); // jk initialize to -1 for "not computed"

	Real const four_pi = 4.0f * Real( numeric::constants::d::pi );

	int num_ones = 0;
	for ( Size ii=1; ii <= pose.total_residue(); ++ii ) {
		Residue const & rsd( pose.residue(ii) );
		rsd_sasa[ ii ] = 0.0;

#ifdef FILE_DEBUG
		TR << "sasa dots: res: " << rsd.name3() << pose.pdb_info()->number(ii) << std::endl;
#endif
		for ( Size iia = 1; iia <= rsd.natoms(); ++iia ) {

			AtomID const id( iia, ii );
			if ( ! atom_subset[ id ] )
				continue; // jk skip this atom if not part of the subset

			Real iia_rad = radii[ rsd.atom(iia).type() ];
			if ( include_probe_radius_in_atom_radii ) iia_rad += probe_radius; // this is the default, not used for sasapack

			//j to get SASA:
			//j - count the number of 1's
			//j - figure fraction that they are
			//j - multiply by 4*pi*r_sqared

			int ctr = 0;
			utility::vector1< ObjexxFCL::ubyte > const & iia_masks( atom_masks[ id ] );
			for ( int bb = 1; bb <= num_bytes; ++bb ) {
				ctr += bit_count[ iia_masks[bb] ]; // atm_masks(bb,iia,ii)
			}
			num_ones += ctr;

			Real const fraction_ones = static_cast< Real >( ctr ) / maskbits; //ronj this is equivalent to l(buried)/l(total)
			Real const total_sa = four_pi * ( iia_rad * iia_rad ); //ronj total atom surface area
			Real const area_exposed = ( 1.0f - fraction_ones ) * total_sa; //ronj ones must indicate buried area if we're subtracting from 1.0

#ifdef FILE_DEBUG
			std::cout << "atom: " << rsd.atom_name( iia ) << ", rad: " << ObjexxFCL::format::F(4,2,radii[ rsd.atom(iia).type() ])
				<< ", covered: " << ObjexxFCL::format::I(3,ctr) << ", counts: "; // use std::cout NOT the TR
			print_dot_bit_string( atom_masks[ id ] );
#endif
			atom_sasa[ id ] = area_exposed;
			// jk Water SASA doesn't count toward the residue's SASA
			if ( ! rsd.atom_type(iia).is_h2o() &&
				 ! rsd.atom_type(iia).is_virtual()) {
				rsd_sasa[ ii ] += area_exposed;
				total_sasa += area_exposed;
			}

		} // iia
#ifdef FILE_DEBUG
		std::cout << "num_ones: " << num_ones << ", sasa: " << rsd_sasa[ ii ] << std::endl; // use std::cout NOT the TR
#endif
		num_ones = 0;

	} // ii

#ifdef FILE_DEBUG
	/*std::cout << std::endl;
	for ( core::Size ii=1; ii <= pose.total_residue(); ++ii ) {
		conformation::Residue const & rsd = pose.residue( ii );
		TR << "residue " << rsd.name3() << pose.pdb_info()->number(ii) << " atom_sasas: [ ";
		for ( Size at=1; at <= rsd.natoms(); ++at ) {
			core::id::AtomID atid( at, ii );
			TR << utility::trim( rsd.atom_name( at ) ) << ":" << atom_sasa[ atid ] << ", ";
		}
		TR << "], total SASA: " << rsd_sasa[ ii ] << std::endl;
	}*/
#endif

	return total_sasa;
}


void
calc_atom_masks(
	core::conformation::Residue const & irsd,
	core::conformation::Residue const & jrsd,
	Real const probe_radius,
	Real const cutoff_distance,
	utility::vector1< Real > const & radii,
	id::AtomID_Map< bool > const & atom_subset,
	core::id::AtomID_Map< utility::vector1< ObjexxFCL::ubyte > > & atom_masks
){

	using core::id::AtomID;
	using core::conformation::Atom;


	Size const ii = irsd.seqpos();
	Size const jj = jrsd.seqpos();

	// use distance rather than distance_squared since the nbr_radii might be negative
	//ronj what residue would have a negative nbr radius?
	//ronj looks like GB_placeholder, SUCK, and VIRT residue types can have negative nbr radii
	Real distance_ij_nbratoms = irsd.atom( irsd.nbr_atom() ).xyz().distance( jrsd.atom( jrsd.nbr_atom() ).xyz() );
	if ( distance_ij_nbratoms > ( irsd.nbr_radius() + jrsd.nbr_radius() + cutoff_distance ) )
		return;

	for ( Size iia=1; iia <= irsd.natoms(); ++iia ) {

		//ronj check to see if this atom is in the subset of atoms we're considering. if not, continue to the next one.
		//ronj have to translate the index 'iia' into an AtomID to check if we're in the subset.
		AtomID const iia_atomid( iia, ii );
		if ( ! atom_subset[ iia_atomid ] )
			continue; // jk skip this atom if not part of the subset

		//ronj convert the atom index into an Atom 'iia_atom' to make getting the xyz() and type() easier
		Atom const & iia_atom( irsd.atom( iia ) );
		Vector const & iia_atom_xyz = iia_atom.xyz();
		Real const iia_atom_radius = radii[ iia_atom.type() ] + probe_radius;

		for ( Size jja = 1; jja <= jrsd.natoms(); ++jja ) {

			//ronj for all atoms in residue 'j', check to make sure that atom is in the subset of atoms we're considering
			AtomID const jji_atomid( jja, jj );
			if ( ! atom_subset[ jji_atomid ] )
				continue; // jk skip this atom if not part of the subset

			Atom const & jja_atom( jrsd.atom( jja ) );
			Vector const & jja_atom_xyz = jja_atom.xyz();
			Real const jja_atom_radius = radii[ jja_atom.type() ] + probe_radius;

			Real const distance_ijxyz( iia_atom_xyz.distance( jja_atom_xyz ) ); // could be faster w/o sqrt, using Jeff Gray's rsq_min stuff
			if ( distance_ijxyz <= iia_atom_radius + jja_atom_radius ) {

				if ( distance_ijxyz <= 0.0 )
					continue;

				// account for atom j overlapping atom i:
				// jk Note: compute the water SASA, but DON'T allow the water to contribute to the burial of non-water atoms
				int degree_of_overlap, aphi, theta, point, masknum;

				if ( ! jrsd.atom_type( jja ).is_h2o() ) {
					get_overlap( iia_atom_radius, jja_atom_radius, distance_ijxyz, degree_of_overlap );
					#ifdef FILE_DEBUG
						//TR << "calculated degree of overlap: " << degree_of_overlap << std::endl;
						//TR << "calculating orientation of " << jrsd.name3() << jj << " atom " << jrsd.atom_name( jja ) << " on "
						//	<< irsd.name3() << ii << " atom " << irsd.atom_name ( iia ) << std::endl;
					#endif

					get_orientation( iia_atom_xyz, jja_atom_xyz, aphi, theta, distance_ijxyz );
					point = angles( aphi, theta );
					masknum = point * 100 + degree_of_overlap;
					#ifdef FILE_DEBUG
						//TR << "calculated masknum " << masknum << std::endl;
					#endif

					//ronj overlap bit values for all atoms should have been init'd to zero before the main for loops
					utility::vector1< ObjexxFCL::ubyte > & iia_bit_values = atom_masks[ AtomID( iia, ii ) ];

					// iterate bb over all 21 bytes or 168 bits (of which we care about 162)
					// bitwise_or the atoms current values with the values from the database/masks array
					#ifdef FILE_DEBUG
						//TR << "starting bit values for atom " << irsd.name3() << ii << "-" << irsd.atom_name( iia ) << ": ";
						//print_dot_bit_string( iia_bit_values );
						//TR << "mask bit values for atom " << irsd.name3() << ii << "-" << irsd.atom_name( iia ) << ": ";
					#endif
					for ( int bb = 1, m = masks.index( bb, masknum ); bb <= num_bytes; ++bb, ++m ) {
						iia_bit_values[ bb ] = ObjexxFCL::bit::bit_or( iia_bit_values[ bb ], masks[ m ] );

						#ifdef FILE_DEBUG
							//int bit;
							//TR << (bb-1) * 8 << ":";
							//for ( int index=7; index >= 0; index-- ) {
							//	bit = ( ( (int)masks[m] >> index ) & 1 );
							//	TR << bit;
							//}
							//TR << " ";
						#endif

					}
					#ifdef FILE_DEBUG
						//TR << std::endl;
						//TR << "overlap bit values for atom " << irsd.name3() << ii << "-" << irsd.atom_name( iia ) << ": ";
						//print_dot_bit_string( iia_bit_values );
					#endif
				}

				// account for i overlapping j:
				// jk Note: compute the water SASA, but DON'T allow the water to contribute to the burial of non-water atoms
				// ronj I don't think this is necessary since we'll eventually perform this calculation when we start
				// ronj iterating over the j atoms
				if ( !irsd.atom_type(iia).is_h2o() ) {
					get_overlap( jja_atom_radius, iia_atom_radius, distance_ijxyz, degree_of_overlap );
					#ifdef FILE_DEBUG
						//TR << "calculated degree of overlap: " << degree_of_overlap << std::endl;
						//TR << "calculating orientation of " << irsd.name3() << ii << " atom " << irsd.atom_name( iia ) << " on "
						//	<< jrsd.name3() << jj << " atom " << jrsd.atom_name ( jja ) << std::endl;
					#endif

					get_orientation( jja_atom_xyz, iia_atom_xyz, aphi, theta, distance_ijxyz );
					point = angles( aphi, theta );
					masknum = point * 100 + degree_of_overlap;
					#ifdef FILE_DEBUG
						//TR << "calculated masknum " << masknum << std::endl;
					#endif

					utility::vector1< ObjexxFCL::ubyte > & jja_bit_values( atom_masks[ AtomID( jja, jj ) ] );

					// iterate bb over all 21 bytes or 168 bits (of which we care about 162)
					// bitwise_or the atoms current values with the values from the database/masks array
					#ifdef FILE_DEBUG
						//TR << "mask bit values for atom " << jrsd.name3() << jj << "-" << jrsd.atom_name( jja ) << ": ";
					#endif
					for ( int bb = 1, m = masks.index(bb,masknum); bb <= num_bytes; ++bb, ++m ) {
						jja_bit_values[ bb ] = ObjexxFCL::bit::bit_or( jja_bit_values[ bb ], masks[ m ] );

						#ifdef FILE_DEBUG
							//int bit;
							//TR << (bb-1) * 8 << ":";
							//for ( int index=7; index >= 0; index-- ) {
							//	bit = ( ( (int)masks[m] >> index ) & 1 );
							//	TR << bit;
							//}
							//TR << " ";
						#endif

					}
					#ifdef FILE_DEBUG
						//TR << std::endl;
						//TR << "final bit values for atom " << jrsd.name3() << jj << "-" << jrsd.atom_name( jja ) << ": ";
						//print_dot_bit_string( jja_bit_values );
					#endif
				}

			} // distance_ijxyz <= iia_atom_radius + jja_atom_radius
			#ifdef FILE_DEBUG
				//TR << "------" << std::endl;
			#endif
		} // jja
	} // jr
}


///
/// @begin get_angles
///
/// @brief
/// Returns the number of bytes the overlap arrays use for tracking SASA.
/// Adding this in so that the values in the SASA database files can be used in SASA-based scores. (ronj)
///
int get_num_bytes() { return num_bytes; }

///
/// @begin get_angles
///
/// @brief
/// Returns const access to the angles FArray, which contains the information in the SASA database file sampling/SASA-angles.dat.
/// Adding this in so that the values in the SASA database files can be used in SASA-based scores. (ronj)
///
ObjexxFCL::FArray2D_int const & get_angles() {
	input_sasa_dats(); // read sasa datafiles; will immediately return if already done
	return angles;
}

///
/// @begin get_masks
///
/// @brief
/// Returns const access to the masks FArray, which contains the information in the SASA database file sampling/SASA-masks.dat.
/// Adding this in so that the values in the SASA database files can be used in SASA-based scores. (ronj)
///
ObjexxFCL::FArray2D_ubyte const & get_masks() {
	input_sasa_dats(); // read sasa datafiles; will immediately return if already done
	return masks;
}

///
/// @begin sasa.cc::calc_per_atom_hydrophobic_sasa
///
/// @brief
/// Uses the method above to calculate total SASA and then only looks at the hydrophobic contribution. Returns the total
/// hydrophobic SASA for the passed in pose.  This method is being used for a protein surface score being developed by ronj.
/// Note: Uses an atom id mask that ignores H's in the pose - only sees and computes the SASA for heavy atoms in the pose.
/// This is done to keep things fast. Only computes the amount of hSASA per residue, not per atom. Doesn't make sense to
/// calculate a per-atom hSASA. (ronj)
///
Real
calc_per_res_hydrophobic_sasa( pose::Pose const & pose,
	utility::vector1< Real > & rsd_sasa, utility::vector1< Real > & rsd_hydrophobic_sasa, Real const probe_radius,
	bool use_naccess_sasa_radii) {

	// an atomID map is needed for the calc_per_atom_sasa method; it stores the actual calculated sasa for every atom
	core::id::AtomID_Map< core::Real > atom_sasa;
	core::pose::initialize_atomid_map( atom_sasa, pose, (core::Real)0.0 ); // initialize to 0.0 for "not computed"

	// clear and init the passed in vector of residue sasa
	rsd_sasa.clear();
	rsd_sasa.resize( pose.total_residue(), 0.0 );

	// create an atom_subset mask such that only the heavy atoms will have their sasa computed (ignore H's to make it faster)
	id::AtomID_Map< bool > atom_subset;

	//ronj calling init_heavy_only leads to uninit'd values in the map
	//ronj what happens is that a vector1 exists for every position in the pose, and when init_heavy_only is called,
	//ronj the first heavy atom number of atoms at the beginning of the vector are set to 'true' or whatever the value
	//ronj type is and then the non-heavy atom indices of the vector are init'd with random values. This behavior occurs
	//ronj even though the vector at each residue position is resized to either the number of atoms or the number of heavy
	//ronj atoms depending on what function was called. It seems like the resize method is not really resizing the vector1
	//ronj at each position, leading to incorrect values being calculated in this method if it's used. Not about to try to
	//ronj debug a vector1 problem; therefore, I'm init'ing the values of the atom_subset atomid_map myself.
	//id::initialize_heavy_only( atom_subset, pose, true ); // this call leads to uninit'd values in the map, causing random results
	atom_subset.clear();
	atom_subset.resize( pose.n_residue() );
	for ( Size ii=1; ii <= pose.n_residue(); ++ii ) {
		atom_subset.resize( ii, pose.residue_type(ii).natoms(), false );
		for ( Size jj = 1; jj <= pose.residue_type(ii).nheavyatoms(); ++jj ) {
			atom_subset[ ii ][ jj ] = true;
		}
	}

	core::Real total_sasa = 0.0;
	total_sasa = core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius, false /* no big polar H */, atom_subset, use_naccess_sasa_radii );
	//#ifdef FILE_DEBUG
	TR.Debug << "total_sasa: " << total_sasa << std::endl;
	//#endif

	// now we have to figure out how much hydrophobic sasa each atom/residue has
	core::Real total_hydrophobic_sasa = 0.0;
	core::Real res_hsasa = 0.0;

	rsd_hydrophobic_sasa.clear();
	rsd_hydrophobic_sasa.resize( pose.total_residue(), 0.0 );

	for ( core::Size ii=1; ii <= pose.total_residue(); ++ii ) {
		conformation::Residue const & rsd = pose.residue( ii );
		res_hsasa = 0.0;

		for ( Size at=1; at <= rsd.nheavyatoms(); ++at ) {
			core::id::AtomID atid( at, ii );

			// exclude hydrogens from consideration of hydrophobic SASA.
			// they should be excluded already because of the atom_subset mask, but just in case.
			if ( rsd.atom_type( at ).is_hydrogen() )
				continue;

			if ( rsd.atom_type( at ).element() == "C" || rsd.atom_type( at ).element() == "S" ) {
				res_hsasa += atom_sasa[ atid ];
			}
		}

		rsd_hydrophobic_sasa[ ii ] = res_hsasa;
		total_hydrophobic_sasa += res_hsasa;
	}

#ifdef FILE_DEBUG
	for ( core::Size ii=1; ii <= pose.total_residue(); ++ii ) {
		conformation::Residue const & rsd = pose.residue( ii );
		TR << "residue " << rsd.name3() << ii << " atom_sasas: [ ";
		for ( Size at=1; at <= rsd.natoms(); ++at ) {
			core::id::AtomID atid( at, ii );
			TR << rsd.atom_type( at ).element() << ":" << atom_sasa[ atid ] << ", ";
		}
		TR << "], total hydrophobic SASA: " << rsd_hydrophobic_sasa[ ii ] << std::endl;
	}
#endif

	return total_hydrophobic_sasa;
}

#ifdef FILE_DEBUG
///
/// @begin sasa.cc::print_dot_bit_string
///
/// @brief
/// helper method I was using to try to confirm that the dots are being overlapped and bits are being set correctly (ronj).
///
void print_dot_bit_string( utility::vector1< ObjexxFCL::ubyte > & values ) {
	for ( int bb = 1; bb <= num_bytes; ++bb ) {
		int bit;
//#ifdef FILE_DEBUG
		if ( (bb-1)*8 % 16 == 0 ) std::cout << (bb-1) * 8 << ":";
//#endif
		for ( int index=7; index >= 0; index-- ) {
			bit = ( ( (int)values[ bb ] >> index ) & 1 );
//#ifdef FILE_DEBUG
			std::cout << bit;
//#endif
		}
//#ifdef FILE_DEBUG
		std::cout << " ";
//#endif
	}
//#ifdef FILE_DEBUG
	std::cout << std::endl;
//#endif
}
#endif
} // namespace scoring
} // namespace core
