// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody_design/LeGrandSasa.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <core/scoring/sasa/LeGrandSasa.hh>
#include <core/scoring/sasa/util.hh>

#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <basic/database/open.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
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

// C++ Headers
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>

static basic::Tracer TR( "core.scoring.sasa.LeGrandSasa" );

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// ObjexxFCL serialization headers
#include <utility/serialization/ObjexxFCL/FArray2D.srlz.hh>
#include <utility/serialization/ObjexxFCL/ubyte.srlz.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace sasa {


using namespace ObjexxFCL::format;
using namespace core;
using utility::vector1;

LeGrandSasa::LeGrandSasa(Real probe_radius, SasaRadii radii_set):
	SasaMethod(probe_radius, radii_set)
{
	init();
}

LeGrandSasa::~LeGrandSasa()= default;

void
LeGrandSasa::init() {
	num_bytes_ = 21;
	num_phi_ = 64;
	num_theta_ = 64;
	num_overlaps_ = 100;
	num_orientations_ = 162;
	maskbits_ = 162;

	angles_ = ObjexxFCL::FArray2D_int(num_phi_, num_theta_);
	masks_ = ObjexxFCL::FArray2D_ubyte(num_bytes_, num_overlaps_ * num_orientations_);

	read_masks();
	read_angles();
}

std::string
LeGrandSasa::get_name() const {
	return "LeGrandSasa";
}

Real
LeGrandSasa::calculate(
	const pose::Pose& pose,
	const id::AtomID_Map<bool> & atom_subset,
	id::AtomID_Map<Real> & atom_sasa,
	vector1<Real> & rsd_sasa
) {
	using core::conformation::Residue;
	using core::conformation::Atom;
	using core::id::AtomID;

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

	if ( pose.size() < 1 ) return 0.0; // nothing to do

	Real const big_polar_H_radius( 1.08 ); // increase radius of polar hydrogens, eg when doing unsatisfied donorH check

	//j setup the radii array, indexed by the atom type int. atom index for looking up an extra data type stored in the AtomTypes
	//ronj reads the values out of the database file sasa_radii.txt in the extras folder of atom_type_sets and stores the values
	//ronj for each atom type into the radii array. each index of the radii array corresponds to some atom type.
	core::chemical::AtomTypeSet const & atom_type_set = pose.residue(1).atom_type_set();

	utility::vector1< Real > radii( atom_type_set.n_atomtypes() );

	for ( core::Size ii=1; ii <= radii.size(); ++ii ) {
		core::chemical::AtomType const & at( atom_type_set[ii] );

		if ( radii_set_ == LJ ) {
			radii[ii] = atom_type_set[ii].lj_radius();
		} else {
			core::Size SASA_RADIUS_INDEX = atom_type_set.extra_parameter_index(get_sasa_radii_parameter_name(radii_set_));
			radii[ii] = atom_type_set[ii].extra_parameter( SASA_RADIUS_INDEX );
		}

		if ( use_big_polar_H_ && at.is_polar_hydrogen() && big_polar_H_radius > radii[ii] ) {
			TR << "Using " << big_polar_H_radius << " instead of " << radii[ii] << " for atomtype " << at.name() << " in sasa calculation!\n";
			radii[ii] = big_polar_H_radius;
		}
	}

	// create an AtomID_Map which will map atom ids to the mask for that atom. right now, just make all atom ids point
	// to a mask of all zeros. later we'll have to set different masks.
	core::id::AtomID_Map< utility::vector1< ObjexxFCL::ubyte > > atom_masks;
	utility::vector1< ObjexxFCL::ubyte > zero_mask( num_bytes_, ObjexxFCL::ubyte( 0 ) );
	core::pose::initialize_atomid_map( atom_masks, pose, zero_mask ); // calls AtomID_Map.Pose.hh::initialize()

	// identify the maxium radius we have for all atoms in the pose.  this will be used to set a cutoff distance for use
	// in skipping residue pairs if their nbr atoms are too far apart.
	core::Real max_radius = 0.0;
	for ( Size ii=1; ii <= pose.size(); ++ii ) {
		Residue const & rsd( pose.residue( ii ) );
		// iterate over all the atoms in residue ii and
		for ( Size jj=1; jj <= rsd.natoms(); ++jj ) {
			max_radius = std::max( max_radius, radii[ rsd.atom(jj).type() ] );
			//ronj assert that the sum of the first 8 bits and the last 8 bits (or first byte and last byte) in the vector
			//ronj of masks for this atom id is less than 0.001.  what is the point of this assertion?
			debug_assert( atom_masks[ AtomID(jj,ii) ][1] + atom_masks[ AtomID(jj,ii) ][num_bytes_] < 1e-3 );
		}
	}

	core::Real cutoff_distance = 2 * ( max_radius + probe_radius_ );

	//j now do calculations: get the atom_masks by looping over all_atoms x all_atoms
	for ( Size ii=1; ii <= pose.size(); ++ii ) {
		Residue const & irsd( pose.residue( ii ) );

		//ronj for the other 'j' residue, only iterate over residues which have indexes > residue 'i'
		for ( Size jj=ii; jj <= pose.size(); ++jj ) {
			Residue const & jrsd( pose.residue( jj ) );

			calc_atom_masks(
				irsd, jrsd,
				probe_radius_, cutoff_distance, radii,
				atom_subset,
				atom_masks );
		}
	}

	//j calculate the residue and atom sasa
	Real total_sasa( 0.0 );

	//JAB Carrying this over for compatability
	atom_sasa.clear();
	core::pose::initialize_atomid_map( atom_sasa, pose, (Real) -1.0 ); // jk initialize to -1 for "not computed"

	Real const four_pi = 4.0f * Real( numeric::constants::d::pi );

	int num_ones = 0;
	for ( Size ii=1; ii <= pose.size(); ++ii ) {
		Residue const & rsd( pose.residue(ii) );

		//TR.Debug << "sasa dots: res: " << rsd.name3() << pose.pdb_info()->number(ii) << std::endl;
		for ( Size iia = 1; iia <= rsd.natoms(); ++iia ) {

			AtomID const id( iia, ii );
			if ( ! atom_subset[ id ] ) {
				continue; // jk skip this atom if not part of the subset
			}

			Real iia_rad = radii[ rsd.atom(iia).type() ];
			if ( include_probe_radius_ ) iia_rad = iia_rad + probe_radius_; // this is the default, not used for sasapack

			//j to get SASA:
			//j - count the number of 1's
			//j - figure fraction that they are
			//j - multiply by 4*pi*r_sqared

			int ctr = 0;
			utility::vector1< ObjexxFCL::ubyte > const & iia_masks( atom_masks[ id ] );
			for ( int bb = 1; bb <= num_bytes_; ++bb ) {
				ctr = ctr + bit_count[ iia_masks[bb] ]; // atm_masks(bb,iia,ii)
			}
			num_ones = num_ones + ctr;

			Real const fraction_ones = static_cast< Real >( ctr ) / maskbits_; //ronj this is equivalent to l(buried)/l(total)
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
					! rsd.atom_type(iia).is_virtual() ) {
				total_sasa  = total_sasa + area_exposed;
				rsd_sasa[ ii ] = rsd_sasa[ii] + area_exposed;
			}

		} // iia
		TR.Debug << "num_ones: " << num_ones << ", sasa: " << std::endl; // use std::cout NOT the TR
		num_ones = 0;

	} // ii

	return total_sasa;
}

void
LeGrandSasa::read_angles() {
	//ronj iterate over all 64 lines of the sampling/SASA-angles.database file and save them to the FArray 'angles'
	//ronj these lines have 64 space-delimited ints which need to be stored
	utility::io::izstream angles_stream( basic::database::full_name( "sampling/SASA-angles.dat" ) );

	for ( int ii = 1; ii <= num_phi_; ++ii ) {
		for ( int jj = 1; jj <= num_theta_; ++jj ) {
			angles_stream >> angles_( ii, jj );
		}
		angles_stream >> skip;
	}
	angles_stream.close();
}

void
LeGrandSasa::read_masks() {
	//j inputting the masks. they are 21 ubytes long, 162x100 (see header). expects file to be complete
	utility::io::izstream masks_stream( basic::database::full_name("sampling/SASA-masks.dat" ) );

	//ronj iterate over all 16200 bit masks in the database file and save them to the FArray 'masks'
	//ronj the lines have 21 space-delimited DECIMAL values which we need to convert to binary values that are 8-bits long
	//ronj so save the int (or short) somewhere and then cast that value into an unsigned char

	for ( int ii=1; ii <= num_overlaps_ * num_orientations_; ++ii ) {
		short tmp; // a 'short' is a int type that is usually 16-bits
		for ( int jj = 1; jj <= num_bytes_; ++jj ) {
			//ronj operator>> of izstream (aka the extraction operator) performs an input operation on a stream generally
			//ronj involving some interpretation of the data
			masks_stream >> tmp;
			masks_(jj,ii) = static_cast< unsigned char >( tmp );
		}
		masks_stream >> skip; // calls skip() on the izstream which skips the rest of the line and the line terminator
	}
	masks_stream.close();
}


void
LeGrandSasa::get_overlap(
	Real const radius_a,
	Real const radius_b,
	Real const distance_ijxyz,
	int & degree_of_overlap
) const {
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
		degree_of_overlap = static_cast< int >( (1.0f - cosine_theta_iq) * 50 ) + 1;

		if ( degree_of_overlap > 100 ) {
			degree_of_overlap = 100;
		} else if ( degree_of_overlap < 0 ) {
			//j We already hopefully accounted for this possibility by requiring that distance_ijxyz > epsilon,
			//j but in case not we don't want a potential bug to go unnoticed.
			TR << "Problem in calculating overlap between atoms. Terminating calculation." << std::endl
				<< "radius_a: " << SS( radius_a ) << ", radius_b: " << SS( radius_b )
				<< ", distance_ijxyz: " << SS( distance_ijxyz ) << ", cosine theta: " << SS( cosine_theta_iq ) << std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}
	}
}


void
LeGrandSasa::get_orientation(
	Vector const & a_xyz,
	Vector const & b_xyz,
	int & phi_index,
	int & theta_index,
	Real distance_ijxyz
) const {
	using namespace numeric::constants::d;
	using numeric::sin_cos_range;

	//pb -- static can cause problems in multi-threading
	Vector ab_diff( ( a_xyz - b_xyz ) / distance_ijxyz );

	//j figuring phi
	//ronj sin_cos_range adjust the range of the passed in argument to a valid [-1,1]. acos(x) returns the arccos of x
	//ronj expressed in radians. (arccos is the inverse operation of cos.) the return value is in the range [0,pi].
	Real phi = std::acos( sin_cos_range( ab_diff(3) ) );

	phi = phi * ( num_phi_ / pi_2 );
	phi_index = static_cast< int >( phi );
	++phi_index; // for fortran goes from 1 to n
	if ( phi_index > num_phi_ ) {
		phi_index = 1; // reset back to index of 1?
	}

	Real theta = std::atan2( ab_diff(2), ab_diff(1) );

	theta = theta * ( num_theta_ / pi_2 );
	theta_index = static_cast< int >( theta );
	++theta_index; // for fortran goes from 1 to n
	if ( theta_index < 1 ) { // as in the case when atan2 returns a negative...
		theta_index  = theta_index + num_theta_; // let t = -32; -32 + 64 = 32
	} else if ( theta_index > num_theta_ ) {
		theta_index = 1; // reset back to index of 1?
	}
}


void
LeGrandSasa::get_2way_orientation(
	Vector const & a_xyz,
	Vector const & b_xyz,
	int & phi_a2b_index,
	int & theta_a2b_index,
	int & phi_b2a_index,
	int & theta_b2a_index,
	Real distance_ijxyz
) const {
	using namespace numeric::constants::d;
	using numeric::sin_cos_range;

	//pb -- static can cause problems in multi-threading
	Vector ab_diff( ( a_xyz - b_xyz ) / distance_ijxyz );

	//j figuring phi, for a2b
	//ronj sin_cos_range adjust the range of the passed in argument to a valid [-1,1]. acos(x) returns the arccos of x
	//ronj expressed in radians. (arccos is the inverse operation of cos.) the return value is in the range [0,pi].
	Real phi_a2b = std::acos( sin_cos_range( ab_diff(3) ) );

	phi_a2b = phi_a2b * ( num_phi_ / pi_2 );
	phi_a2b_index = static_cast< int >( phi_a2b );
	++phi_a2b_index; // for fortran goes from 1 to n
	if ( phi_a2b_index > num_phi_ ) {
		phi_a2b_index = 1; // reset back to index of 1?
	}

	// can take a shortcut to get phi_b2a
	Real phi_b2a = num_phi_ / 2 - phi_a2b;
#ifdef FILE_DEBUG
	Real temp_phi_b2a = phi_b2a;
#endif
	phi_b2a_index = static_cast< int >( phi_b2a );
	++phi_b2a_index; // for fortran goes from 1 to n
	if ( phi_b2a_index > num_phi_ ) {
		phi_b2a_index = 1;
	}

	Real theta_a2b = std::atan2( ab_diff(2), ab_diff(1) );

	theta_a2b = theta_a2b * ( num_theta_ / pi_2 );
	theta_a2b_index = static_cast< int >( theta_a2b );
	++theta_a2b_index; // for fortran goes from 1 to n
	if ( theta_a2b_index < 1 ) { // as in the case when atan2 returns a negative...
		theta_a2b_index = (core::Size)(theta_a2b + num_theta_); // let t = -32; -32 + 64 = 32
	} else if ( theta_a2b_index > num_theta_ ) {
		theta_a2b_index = 1; // reset back to index of 1?
	}

	// can use a shortcut to get theta_b2a
	Real theta_b2a = num_theta_ / 2.0f + theta_a2b;
	if ( theta_b2a > num_theta_ / 2.0f ) {
		theta_b2a -= num_theta_;
	}

	theta_b2a_index = static_cast< int >( theta_b2a );
	++theta_b2a_index; // for fortran goes from 1 to n
	if ( theta_b2a_index < 1 ) {
		theta_b2a_index = theta_b2a_index + num_theta_;
	} else if ( theta_b2a_index > num_theta_ ) {
		theta_b2a_index = 1;
	}
}


void
LeGrandSasa::calc_atom_masks(
	core::conformation::Residue const & irsd,
	core::conformation::Residue const & jrsd,
	Real const probe_radius,
	Real const cutoff_distance,
	utility::vector1< Real > const & radii,
	id::AtomID_Map< bool > const & atom_subset,
	core::id::AtomID_Map< utility::vector1< ObjexxFCL::ubyte > > & atom_masks
) const {
	using core::id::AtomID;
	using core::conformation::Atom;

	Size const ii = irsd.seqpos();
	Size const jj = jrsd.seqpos();

	// use distance rather than distance_squared since the nbr_radii might be negative
	//ronj what residue would have a negative nbr radius?
	//ronj looks like GB_placeholder, SUCK, and VIRT residue types can have negative nbr radii
	Real distance_ij_nbratoms = irsd.atom( irsd.nbr_atom() ).xyz().distance( jrsd.atom( jrsd.nbr_atom() ).xyz() );
	if ( distance_ij_nbratoms > ( irsd.nbr_radius() + jrsd.nbr_radius() + cutoff_distance ) ) {
		return;
	}

	for ( Size iia=1; iia <= irsd.natoms(); ++iia ) {

		//ronj check to see if this atom is in the subset of atoms we're considering. if not, continue to the next one.
		//ronj have to translate the index 'iia' into an AtomID to check if we're in the subset.
		AtomID const iia_atomid( iia, ii );
		if ( ! atom_subset[ iia_atomid ] ) continue; // jk skip this atom if not part of the subset

		//ronj convert the atom index into an Atom 'iia_atom' to make getting the xyz() and type() easier
		Atom const & iia_atom( irsd.atom( iia ) );
		Vector const & iia_atom_xyz = iia_atom.xyz();
		Real const iia_atom_radius = radii[ iia_atom.type() ] + probe_radius;

		for ( Size jja = 1; jja <= jrsd.natoms(); ++jja ) {

			//ronj for all atoms in residue 'j', check to make sure that atom is in the subset of atoms we're considering
			AtomID const jji_atomid( jja, jj );
			if ( ! atom_subset[ jji_atomid ] ) continue; // jk skip this atom if not part of the subset

			Atom const & jja_atom( jrsd.atom( jja ) );
			Vector const & jja_atom_xyz = jja_atom.xyz();
			Real const jja_atom_radius = radii[ jja_atom.type() ] + probe_radius;

			Real const distance_ijxyz( iia_atom_xyz.distance( jja_atom_xyz ) ); // could be faster w/o sqrt, using Jeff Gray's rsq_min stuff
			if ( distance_ijxyz > iia_atom_radius + jja_atom_radius ) {
				TR.Debug << "------" << std::endl;
				continue;
			}

			if ( distance_ijxyz <= 0.0 ) continue;

			// account for atom j overlapping atom i:
			// jk Note: compute the water SASA, but DON'T allow the water to contribute to the burial of non-water atoms
			int degree_of_overlap, aphi, theta, point, masknum;

			if ( ! jrsd.atom_type( jja ).is_h2o() ) {
				get_overlap( iia_atom_radius, jja_atom_radius, distance_ijxyz, degree_of_overlap );

				TR.Debug << "calculated degree of overlap: " << degree_of_overlap << std::endl
					<< "calculating orientation of " << jrsd.name3() << jj << " atom " << jrsd.atom_name( jja ) << " on "
					<< irsd.name3() << ii << " atom " << irsd.atom_name ( iia ) << std::endl;

				get_orientation( iia_atom_xyz, jja_atom_xyz, aphi, theta, distance_ijxyz );
				point = angles_( aphi, theta );
				masknum = point * 100 + degree_of_overlap;

				TR.Debug << "calculated masknum " << masknum << std::endl;

				//ronj overlap bit values for all atoms should have been init'd to zero before the main for loops
				utility::vector1< ObjexxFCL::ubyte > & iia_bit_values = atom_masks[ AtomID( iia, ii ) ];

				// iterate bb over all 21 bytes or 168 bits (of which we care about 162)
				// bitwise_or the atoms current values with the values from the database/masks array
#ifdef FILE_DEBUG
				//TR << "starting bit values for atom " << irsd.name3() << ii << "-" << irsd.atom_name( iia ) << ": ";
				//print_dot_bit_string( iia_bit_values );
				//TR << "mask bit values for atom " << irsd.name3() << ii << "-" << irsd.atom_name( iia ) << ": ";
#endif
				for ( int bb = 1, m = masks_.index( bb, masknum ); bb <= num_bytes_; ++bb, ++m ) {
					iia_bit_values[ bb ] = ObjexxFCL::bit::bit_or( iia_bit_values[ bb ], masks_[ m ] );

#ifdef FILE_DEBUG
					//int bit;
					//TR << (bb-1) * 8 << ":";
					//for ( int index=7; index >= 0; index-- ) {
					//	bit = ( ( (int)masks_[m] >> index ) & 1 );
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
			if ( !irsd.atom_type(iia).is_h2o() ) {
				get_overlap( jja_atom_radius, iia_atom_radius, distance_ijxyz, degree_of_overlap );
				TR.Debug << "calculated degree of overlap: " << degree_of_overlap << std::endl
					<< "calculating orientation of " << irsd.name3() << ii << " atom " << irsd.atom_name( iia ) << " on "
					<< jrsd.name3() << jj << " atom " << jrsd.atom_name ( jja ) << std::endl;

				get_orientation( jja_atom_xyz, iia_atom_xyz, aphi, theta, distance_ijxyz );
				point = angles_( aphi, theta );
				masknum = point * 100 + degree_of_overlap;

				TR.Debug << "calculated masknum " << masknum << std::endl;

				utility::vector1< ObjexxFCL::ubyte > & jja_bit_values( atom_masks[ AtomID( jja, jj ) ] );

				// iterate bb over all 21 bytes or 168 bits (of which we care about 162)
				// bitwise_or the atoms current values with the values from the database/masks array
				TR.Debug << "mask bit values for atom " << jrsd.name3() << jj << "-" << jrsd.atom_name( jja ) << ": ";

				for ( int bb = 1, m = masks_.index(bb,masknum); bb <= num_bytes_; ++bb, ++m ) {
					jja_bit_values[ bb ] = ObjexxFCL::bit::bit_or( jja_bit_values[ bb ], masks_[ m ] );

#ifdef FILE_DEBUG
					//int bit;
					//TR << (bb-1) * 8 << ":";
					//for ( int index=7; index >= 0; index-- ) {
					//	bit = ( ( (int)masks_[m] >> index ) & 1 );
					//	TR << bit;
					//}
					//TR << " ";
#endif
				}
#ifdef FILE_DEBUG
				//TR.Debug << std::endl
				//TR << "final bit values for atom " << jrsd.name3() << jj << "-" << jrsd.atom_name( jja ) << ": ";
				//print_dot_bit_string( jja_bit_values );
#endif
			}
			TR.Debug << "------" << std::endl;
		}
	}
}


ObjexxFCL::FArray2D_int const &
LeGrandSasa::get_angles() const {
	return angles_;
}


ObjexxFCL::FArray2D_ubyte const &
LeGrandSasa::get_masks() const {
	return masks_;
}


/// @brief
/// helper method I was using to try to confirm that the dots are being overlapped and bits are being set correctly (ronj).
///
void
LeGrandSasa::print_dot_bit_string( utility::vector1< ObjexxFCL::ubyte > & values ) const{
	for ( int bb = 1; bb <= num_bytes_; ++bb ) {
		if ( (bb-1)*8 % 16 == 0 ) std::cout << (bb-1) * 8 << ":";
		for ( int index=7; index >= 0; index-- ) {
			int bit = ( ( (int)values[ bb ] >> index ) & 1 );
			TR << bit;
		}
		TR << " ";
	}
	TR << std::endl;
}

} //sasa
} //scoring
} //core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::sasa::LeGrandSasa::LeGrandSasa() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::sasa::LeGrandSasa::save( Archive & arc ) const {
	arc( cereal::base_class< core::scoring::sasa::SasaMethod >( this ) );
	arc( CEREAL_NVP( num_bytes_ ) ); // int
	arc( CEREAL_NVP( num_phi_ ) ); // int
	arc( CEREAL_NVP( num_theta_ ) ); // int
	arc( CEREAL_NVP( num_overlaps_ ) ); // int
	arc( CEREAL_NVP( num_orientations_ ) ); // int
	arc( CEREAL_NVP( maskbits_ ) ); // int
	arc( CEREAL_NVP( angles_ ) ); // ObjexxFCL::FArray2D<int>
	arc( CEREAL_NVP( masks_ ) ); // ObjexxFCL::FArray2D<ObjexxFCL::ubyte>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::sasa::LeGrandSasa::load( Archive & arc ) {
	arc( cereal::base_class< core::scoring::sasa::SasaMethod >( this ) );
	arc( num_bytes_ ); // int
	arc( num_phi_ ); // int
	arc( num_theta_ ); // int
	arc( num_overlaps_ ); // int
	arc( num_orientations_ ); // int
	arc( maskbits_ ); // int
	arc( angles_ ); // ObjexxFCL::FArray2D<int>
	arc( masks_ ); // ObjexxFCL::FArray2D<ObjexxFCL::ubyte>
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::sasa::LeGrandSasa );
CEREAL_REGISTER_TYPE( core::scoring::sasa::LeGrandSasa )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_sasa_LeGrandSasa )
#endif // SERIALIZATION
