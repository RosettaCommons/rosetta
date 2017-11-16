// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody_design/util.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <basic/database/open.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/sasa/SasaCalc.hh>
#include <core/scoring/sasa/util.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

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
#include <core/pose/util.hh>

// C++ Headers
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include<algorithm>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>

static basic::Tracer TR( "core.scoring.sasa" );

using namespace ObjexxFCL::format;

namespace core {
namespace scoring {
namespace sasa {
using namespace core;
using namespace ObjexxFCL::format;

///////////////////////////////////////////////////////////////////////////////////////////////
///  Convenience Functions
///
///

/// @brief Absolute per residue sidechain SASA
/// @details Added by JKLeman (julia.koehler1982@gmail.com)
///   GXG tripeptide values for sidechain SASA are taken from
///   http://www.proteinsandproteomics.org/content/free/tables_1/table08.pdf
utility::vector1< Real > per_res_sc_sasa( const pose::Pose & pose ) {

	// get per-residue SASA
	SasaCalc calc = SasaCalc();
	calc.calculate( pose );

	// get residue SASA
	return calc.get_residue_sasa_sc();

} // per residue sasa

/// @brief Relative per residue sidechain SASA
/// @details Added by JKLeman (julia.koehler1982@gmail.com)
///   GXG tripeptide values for sidechain SASA are taken from
///   http://www.proteinsandproteomics.org/content/free/tables_1/table08.pdf
utility::vector1< Real > rel_per_res_sc_sasa( const pose::Pose & pose ) {

	using namespace core::pose;

	// get residue SASA
	utility::vector1< Real > res_sasa = per_res_sc_sasa( pose );

	// create map of sidechain SASAs
	// http://www.proteinsandproteomics.org/content/free/tables_1/table08.pdf
	std::map< char, Real > sc_sasa;

	sc_sasa[ 'A' ] = 67;
	sc_sasa[ 'R' ] = 196;
	sc_sasa[ 'N' ] = 113;
	sc_sasa[ 'D' ] = 106;
	sc_sasa[ 'C' ] = 104;
	sc_sasa[ 'Q' ] = 144;
	sc_sasa[ 'E' ] = 138;
	sc_sasa[ 'G' ] = 1;
	sc_sasa[ 'H' ] = 151;
	sc_sasa[ 'I' ] = 140;
	sc_sasa[ 'L' ] = 137;
	sc_sasa[ 'K' ] = 167;
	sc_sasa[ 'M' ] = 160;
	sc_sasa[ 'F' ] = 175;
	sc_sasa[ 'P' ] = 105;
	sc_sasa[ 'S' ] = 80;
	sc_sasa[ 'T' ] = 102;
	sc_sasa[ 'W' ] = 217;
	sc_sasa[ 'Y' ] = 187;
	sc_sasa[ 'V' ] = 117;

	// create vector with relative sasa
	utility::vector1< Real > rel_sasa;
	for ( Size i = 1; i <= nres_protein( pose ); ++i ) {
		// debug
		//  if ( pose.residue( i ).name3() == "GLY" ) {
		//   TR << "res_sasa " << i << " " << res_sasa[ i ] << std::endl;
		//  }
		rel_sasa.push_back( res_sasa[ i ] / sc_sasa[ pose.residue( i ).name1() ] );
	}

	return rel_sasa;
} // relative per residue sidechain sasa

/// @brief Is residue exposed?
/// @details Added by JKLeman (julia.koehler1982@gmail.com)
///   Uses the function rel_per_res_sc_sasa above
///   THIS IS EXPENSIVE, BE AWARE!!! IF YOU NEED TO RUN IT OVER THE ENTIRE
///   PROTEIN, USE THE rel_per_res_sc_sasa FUNCTION INSTEAD!
bool is_res_exposed( const pose::Pose & pose, core::Size resnum ) {

	// get relative sidechain SASA
	utility::vector1< Real > rel_sasa = rel_per_res_sc_sasa( pose );

	// if relative sc SASA > 0.5, residue is exposed
	// otherwise residue is buried
	return ( rel_sasa[ resnum ] >= 0.5 );
} // is residue exposed?


/// @brief Calculate the sidechain and backbone sasa from atom sasa
std::pair<Real, Real>
get_sc_bb_sasa(const pose::Pose & pose, const id::AtomID_Map<Real> & atom_sasa) {

	Real sc_sasa = 0.0;
	Real bb_sasa = 0.0;

	for ( Size i = 1; i <= atom_sasa.size(); ++i ) {
		utility::vector1< Size> bb_atoms = pose.residue_type(i).all_bb_atoms();
		utility::vector1< Size> sc_atoms = pose.residue_type(i).all_sc_atoms();

		//BB Calculation
		for ( Size x = 1; x <= bb_atoms.size(); ++x ) {
			core::id::AtomID atomid(bb_atoms[x], i);
			if ( atom_sasa[atomid] < 0 ) continue; //Non-computed sasa
			bb_sasa += atom_sasa[atomid];
		}

		//SC Calculation
		for ( Size x = 1; x <= sc_atoms.size(); ++x ) {
			core::id::AtomID atomid(sc_atoms[x], i);
			if ( atom_sasa[atomid] < 0 ) continue; // Non-computed sasa
			sc_sasa += atom_sasa[atomid];

		}
	}

	std::pair<Real, Real> split_sasa = std::make_pair(sc_sasa, bb_sasa);
	return split_sasa;
}

std::pair<utility::vector1<Real>, utility::vector1<Real> >
get_sc_bb_sasa_per_res(const pose::Pose & pose, const id::AtomID_Map<Real> & atom_sasa) {

	utility::vector1< Real > sc_sasa(pose.size(), 0.0);
	utility::vector1< Real > bb_sasa(pose.size(), 0.0);

	for ( Size i = 1; i <= atom_sasa.size(); ++i ) {
		utility::vector1< Size> bb_atoms = pose.residue_type(i).all_bb_atoms();
		utility::vector1< Size> sc_atoms = pose.residue_type(i).all_sc_atoms();

		//BB Calculation
		for ( Size x = 1; x <= bb_atoms.size(); ++x ) {
			core::id::AtomID atomid(bb_atoms[x], i);
			if ( atom_sasa[atomid] < 0 ) continue; //Non-computed sasa
			bb_sasa[i] += atom_sasa[atomid];
		}

		//SC Calculation
		for ( Size x = 1; x <= sc_atoms.size(); ++x ) {
			core::id::AtomID atomid(sc_atoms[x], i);
			if ( atom_sasa[atomid] < 0 ) continue; // Non-computed sasa
			sc_sasa[i] += atom_sasa[atomid];

		}
	}

	std::pair<utility::vector1<Real>, utility::vector1<Real> > split_sasa = std::make_pair(sc_sasa, bb_sasa);
	return split_sasa;
}


///////////////////////////////////////////////////////////////////////////////////////////////
///  Enum Management.  Can go in separate file if it gets large.
///
///
///

SasaMethodEnum
get_sasa_method_from_string(std::string method) {
	std::map<std::string, SasaMethodEnum> methods;

	methods["LeGrand"] = LeGrand;

	return methods[method];
}

std::string
get_sasa_radii_parameter_name(SasaRadii radii_set) {
	utility::vector1<std::string> radii_names;

	radii_names.resize(SasaRadii_total);

	radii_names[LJ] = "LJ_RADII";
	radii_names[legacy] = "SASA_RADIUS_LEGACY";
	radii_names[naccess] = "NACCESS_SASA_RADIUS";
	radii_names[reduce] = "REDUCE_SASA_RADIUS";
	radii_names[chothia] = radii_names[naccess];

	return radii_names[radii_set];
}

SasaRadii
get_sasa_radii_set_from_string(std::string radii_set){
	std::map<std::string, SasaRadii> radii_enums;
	radii_enums["chothia"] = chothia;
	radii_enums["naccess"] = naccess;
	radii_enums["reduce"] = reduce;
	radii_enums["legacy"] = legacy;
	radii_enums["LJ"] = LJ;

	return radii_enums[radii_set];
}



///////////////////////////////////////////////////////////////////////////////////////////
//// Functions used for molecular surface approximation
//// LeGrand S, Merz KM. Rapid approximation to molecular surface area via the use of Boolean logic and look-up tables.
////   J Comput Chem 1993;14:349–352.
////
////  Implementation: Jerry Tsai
////  C++ Translation: Jeff Gray
////


// These are needed currently by get_legrand_sasa_masks and get_legrand_sasa_angles:
int const num_bytes = 21;
int const num_phi = 64;
int const num_theta = 64;
int const num_overlaps = 100;
int const num_orientations = 162;

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
/// @brief
/// Returns const access to the masks FArray, which contains the information in the SASA database file sampling/SASA-masks.dat.
/// Adding this in so that the values in the SASA database files can be used in SASA-based scores. (ronj)
///
ObjexxFCL::FArray2D_ubyte const & get_legrand_sasa_masks() {
	input_legrand_sasa_dats(); // read sasa datafiles; will immediately return if already done
	return masks;
}

///
/// @brief
/// Returns const access to the angles FArray, which contains the information in the SASA database file sampling/SASA-angles.dat.
/// Adding this in so that the values in the SASA database files can be used in SASA-based scores. (ronj)
///
ObjexxFCL::FArray2D_int const & get_legrand_sasa_angles() {
	input_legrand_sasa_dats(); // read sasa datafiles; will immediately return if already done
	return angles;
}


/// @brief
/// Reads in the SASA database files sampling/SASA-angles.dat and sampling/SASA-masks.dat into FArrays above.
///
void input_legrand_sasa_dats() {

	// this booleans checks whether we've done this already or not. if we have, then return immediately.
	static bool init = false;
	if ( init ) return;

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
/// @details
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
get_legrand_atomic_overlap( Real const radius_a, Real const radius_b, Real const distance_ijxyz, int & degree_of_overlap ) {

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
/// @brief
/// Gets the orientation of a to b (i to j, see below). Does this by calculating two angles, aphi and theta. (j)
///
/// @details
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
void get_legrand_orientation( Vector const & a_xyz, Vector const & b_xyz, int & phi_index, int & theta_index, Real distance_ijxyz ) {

	using namespace numeric::constants::d;
	using numeric::sin_cos_range;

	//pb -- static can cause problems in multi-threading
	Vector ab_diff( ( a_xyz - b_xyz ) / distance_ijxyz );

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
	if ( phi_index > num_phi ) {
		phi_index = 1; // reset back to index of 1?
	}

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
/// @brief
/// Gets the orientation of a to b (i to j, see below). Does this by calculating two angles, aphi and theta. (j)
///
/// @details
/// ronj This function is the same as the function above but get the orientation of a to b simultaneously with the
/// ronj orientation of b to a.  The same result could be achieved by making two separate get_legrand_2way_orientation() calls
/// ronj but this method does it more efficiently by avoiding an atan2 and acos call.  Instead, once you compute the
/// ronj phi and theta for a on b, you can add/subtrate pi factors to get the phi and theta for b on a.
/// ronj Still not sure how this method returns the correct values, though.
///
void
get_legrand_2way_orientation( Vector const & a_xyz, Vector const & b_xyz, int & phi_a2b_index, int & theta_a2b_index, int & phi_b2a_index, int & theta_b2a_index, Real distance_ijxyz ) {

	using namespace numeric::constants::d;
	using numeric::sin_cos_range;

	//pb -- static can cause problems in multi-threading
	Vector ab_diff( ( a_xyz - b_xyz ) / distance_ijxyz );

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
	if ( phi_a2b_index > num_phi ) {
		phi_a2b_index = 1; // reset back to index of 1?
	}

	// can take a shortcut to get phi_b2a
	Real phi_b2a = num_phi / 2 - phi_a2b;
#ifdef FILE_DEBUG
	Real temp_phi_b2a = phi_b2a;
#endif
	phi_b2a_index = static_cast< int >( phi_b2a );
	++phi_b2a_index; // for fortran goes from 1 to n
	if ( phi_b2a_index > num_phi ) {
		phi_b2a_index = 1;
	}

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
	if ( theta_b2a > num_theta / 2.0f ) {
		theta_b2a -= num_theta;
	}
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



}
}
}
