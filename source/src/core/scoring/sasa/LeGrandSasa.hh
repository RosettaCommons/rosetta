// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file core/scoring/sasa/LeGrandSasa.hh
/// @brief Sasa Method LeGrand 1993
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_scoring_sasa_LEGRANDSASA_HH
#define INCLUDED_core_scoring_sasa_LEGRANDSASA_HH


#include <core/scoring/sasa/SasaMethod.hh>
#include <core/scoring/sasa/LeGrandSasa.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>

// Utility headers


#include <utility/vector1.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/ubyte.hh>

namespace core {
namespace scoring {
namespace sasa {
	using namespace core;
	using utility::vector1;
	using ObjexxFCL::ubyte;
	
///@brief LeGrand SASA approximation method
///  Used by SasaCalc but can be used by itself.
///
///@details LeGrand S, Merz KM. Rapid approximation to molecular surface area via the use of Boolean logic and look-up tables.
///   J Comput Chem 1993;14:349 –352.
///
///  Fortran Implementation: Jerry Tsai
///  C++ Translation: Jeff Gray
///  Cleanup/Bugfixes/OOP: Jared Adolf-Bryfogle
///
class LeGrandSasa : public SasaMethod {
	

public:
	
	LeGrandSasa(Real probe_radius, SasaRadii radii_set);
	virtual ~ LeGrandSasa();
	
	  
	//Real
	//calculate(const pose::Pose & pose);
	
	///@brief Calculate Sasa.  Atoms not calculated have -1 sasa.  This is carried over for compatability purposes.  	
	virtual Real
	calculate(
		const pose::Pose & pose,
		const id::AtomID_Map<bool> & atom_subset,
		id::AtomID_Map<Real> & atom_sasa,
		vector1< Real > & rsd_sasa);
	
	virtual std::string
	get_name() const ;

public:
	
	///@brief
	/// Returns const access to the angles FArray, which contains the information in the SASA database file sampling/SASA-angles.dat.
	/// Adding this in so that the values in the SASA database files can be used in SASA-based scores. (ronj)
	///
	ObjexxFCL::FArray2D_int const &
	get_angles() const;
	
	///@brief
	/// Returns const access to the masks FArray, which contains the information in the SASA database file sampling/SASA-masks.dat.
	/// Adding this in so that the values in the SASA database files can be used in SASA-based scores. (ronj)
	ObjexxFCL::FArray2D_ubyte const &
	get_masks() const;
	
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
	get_overlap( Real const radius_a, Real const radius_b, Real const distance_ijxyz, int & degree_of_overlap ) const;
	
	/// @brief
	/// Gets the orientation of a to b (i to j, see below). Does this by calculating two angles, aphi and theta. (j)
	///
	/// @detailed
	/// This function is used to get two indexes (phi and theta) which are used to get the index of a dot on the
	/// surface of the 'a' sphere. When calculating how much surface area sphere b covers on a, we can get the degree
	/// of overlap from the function above, but it's not necessarily the case that the vector that connects the center
	/// of atom 'a' and atom 'b' goes through one of the predetermined dot locations on the surface of 'a'. In fact,
	/// it's very unlikely that the vector goes through a predetermined dot.  Instead, what is done is the actual point
	/// of intersection (the outermost point of a on the line from the center of 'a' to center of 'b') is converted
	/// to spherical polar coordinates. Then, the values are used to find the location of the closest predetermined
	/// point on the surface of 'a' using a lookup table. So what this function needs to do is convert the
	/// cartesian coordinate of the actual point of intersection into polar coordinates.
	///
	/// To get the spherical, polar coordinates of a cartesian point x,y,z, use these equations:
	/// r = sqrt( x^2 + y^2 + z^2 )
	/// theta = arccos( z / r )
	/// phi = arctan( y / x )
	///
	/// Then, once we have the true phi and theta, we need to translate this into an index (or offset) for the correct
	/// value in the database file. There are 64 phi angle bin and 64 theta bins in the database file sampling/SASA-angles.dat.
	/// We need to convert the phi and theta into indexes for this file by multiplying them by num_phi / 2*pi.
	//  pi_2 is 2*pi
	///
	/// Note: I think phi and theta have been reversed in the function below. The code below uses the following:
	/// phi = arccos( z )
	/// theta = arctan( y / x )
	///
	/// After a couple of weeks trying to write tests for this function, I have been unsuccessful in figuring out why
	/// it's doing what it does.  Despite using the wrong equations, it seems to work.  Comparing the total residue
	/// SASA values calculated by calc_per_atom_sasa() below results in a correlation of 0.98 against what the program
	/// NACCESS finds for the same residues. This test was done on a small 110aa protein.  I also looked at the per-atom
	/// total SASA and the correlation for all atoms (mini v. NACCESS) was approximately 0.94. I'm using exactly the same
	/// van der Waals radii for both programs so I feel like the correlations should be 1.0. Explanations for the
	/// differences can be 1) this method is doing something wrong in calculating the closest surface point, 2) this
	/// method is correct but the masks that are in the database are not aligned to the surface points correctly, 3) the
	/// differences are solely due to the different way that the two program calculate surface area.
	/// (ronj)
	///
	void
	get_orientation( Vector const & a_xyz, Vector const & b_xyz, int & phi_index, int & theta_index, Real distance_ijxyz ) const;
	
	/// @brief
	/// Gets the orientation of a to b (i to j, see below). Does this by calculating two angles, aphi and theta. (j)
	///
	/// @detailed
	/// This function is the same as the function above but get the orientation of a to b simultaneously with the
	/// orientation of b to a.  The same result could be achieved by making two separate get_2way_orientation() calls
	/// but this method does it more efficiently by avoiding an atan2 and acos call.  Instead, once you compute the
	/// phi and theta for a on b, you can add/subtrate pi factors to get the phi and theta for b on a.
	/// Still not sure how this method returns the correct values, though.
	/// (ronj)
	///
	void
	get_2way_orientation( Vector const & a_xyz, Vector const & b_xyz,
		int & phi_a2b_index, int & theta_a2b_index, int & phi_b2a_index, int & theta_b2a_index, Real distance_ijxyz ) const;



	void
	calc_atom_masks(
		conformation::Residue const & irsd,
		conformation::Residue const & jrsd,
		Real const probe_radius,
		Real const cutoff_distance,
		utility::vector1< Real > const & radii,
		id::AtomID_Map< bool > const & atom_subset,
		id::AtomID_Map< utility::vector1< ObjexxFCL::ubyte > > & atom_mask ) const;


private:
	
	///@brief Initialize the class - allows alternate constructors, copy constructors, etc.
	void
	init();
	
	///@brief Read angles db file into private FArray
	void
	read_angles();
	
	///@brief Read masks db file into private FArray
	void
	read_masks();
	
	
	///@brief
	/// helper method to try to confirm that the dots are being overlapped and bits are being set correctly (ronj).
	///
	void
	print_dot_bit_string( utility::vector1< ObjexxFCL::ubyte > & values ) const;
	
	
private:
	
	
	int num_bytes_;
	int num_phi_;
	int num_theta_;
	int num_overlaps_;
	int num_orientations_;
	int maskbits_;	
	
	ObjexxFCL::FArray2D<int> angles_;
	ObjexxFCL::FArray2D<ubyte> masks_;
	
	//std::string angles_db_file_;
	//std::string masks_db_file_;
};


}
}
}


#endif	//#ifndef INCLUDED_protocols/antibody_design_LEGRANDSASA_HH

