// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/ExternalGeomSampler.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_toolbox_match_enzdes_util_ExternalGeomSampler_hh
#define INCLUDED_protocols_toolbox_match_enzdes_util_ExternalGeomSampler_hh

// Unit headers
#include <protocols/toolbox/match_enzdes_util/ExternalGeomSampler.fwd.hh>

// Project headers
#include <core/types.hh>

// Numeric headers
#include <numeric/HomogeneousTransform.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <list>

#include <utility/vector1.hh>


namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

/// @brief The external geom sampler is a class that manages the data necessary
/// to construct the coordinates of the three atoms in the downstream partner.
///
/// @details
/// The external geom sampler holds the data for the six degrees of freedom
/// that describe the downstream atoms relative to the upstream atoms.
/// The naming of the 6 atoms involved and the 6 degrees of freedom involved
/// is as follows:
///
/// Atoms U1, U2, and U3 are on the upstream partner.
/// Atoms D1, D2, and D3 are on the downstream partner.
/// The numbers indicate the distance the atoms have from the "break" between the upstream
/// and downstream partners.
///
/// DOF 1: "tor_U3D1" is the torsion between atoms U3, U2, U1 and D1.
/// DOF 2: "ang_U2D1" is angle between atoms U2, U1 and D1.
/// DOF 3: "dis_U1D1" is the distance between atoms U1 and D1.
/// DOF 4: "tor_U2D2" is the torsion between atoms U2, U1, D1 and D2.
/// DOF 5: "ang_U1D2" is the angle between atoms U1, D1 and D3
/// DOF 6: "tor_U1D3"  is the torsion between atoms U1, D1, D2 and D3
///
/// This naming convention describes a parameter by three letters followed by two atom names.
/// The letter is "d" for distance, "a" for angle, and "t" for torsion.  The two numbers
/// describe the end points of that parameter assuming a continuous span of numbers
/// between.  E.g. t25 describes the geometry of 4 atoms starting at 2 and going to 5,
/// hitting atoms 3 and 4 along the way.
///
/// ASCII art:
///  Fig 1. The 6 atoms that define the geometry between the upstream partner
///  and the downstream partner.  dis_U1D1 is the distance between atoms U1 and D1
///  U3
///     \.
///      U2 --- U1
///
///                \ dis_U1D1
///
///                  D1 --- D2
///                           \.
///                             D3
///
///  Fig 2. The two angles
///  U3
///     \.
///      U2 --- U1
///        ang_U2D1
///                \.
///                   ang_U1D2
///                  D1 --- D2
///                           \.
///                             D3
///  Fig 3. The three torsions
///  U3
///     \  tor_U3D1                tor_U3D1: torsion about U2-U1 axis starting at atom U3 and ending at D1
///      U2 --- U1
///
///                \ tor_U2D2      tor_U2D2: torsion about U1-D1 axis starting at atom U2 and ending at D2
///
///                  D1 --- D2
///                tor_U1D3   \     tor_U1D3: torsion about D1-D2 axis starting at atom U1 and ending at D3
///                             D3
///
/// There are three other parameters needed to construct the location
/// of atoms 4, 5 and 6 from atoms 1, 2 and 3, but these are not
/// considered DOFs as they are not allowed to vary.
/// They are:
/// Extra 1: the distance between atom D1 and atom D2
/// Extra 2: the distance between atom D2 and atom D3
/// Extra 3: the angle between atoms D1, D2 and D3.
///
/// These three parameters must be communicated to the ExternGeomSampler
/// through its functions set_dis_D1D2, set_dis_D2D3 and set_ang_D1D2D3.
///
/// The ExternalGeomSampler holds the set of samples for each of the
/// 6 DOFs.  For each sample, it precomputes a HomogeneousTransform
/// so that during the enumeration of all sample combinations the
/// coordinate computation requires matrix multiplication only: the
/// more expensive transcendental function evaluations evaluated before
/// coordinate construction begins, their results stored, and are thereby
/// avoided in the enumeration loops.
///
/// To build the coordinates of atoms D1, D2 and D3, given the sample id's
/// for the 6 DOFs as the 6-tuple [ i, j, k, l, m, n ]:
/// start with a coordinate frame at atom U1 s.t. the z axis lies on the
/// vector from atom U2 to atom U1, the y axis lies in the plane defined
/// atoms U1, U2 and U3, and the x axis is the cross-product of y and z.
/// Call the start frame "frame1."
/// frame2 is the product: frame1 * transform( HT_tor_U3D1, i );
/// frame3 is the product: frame2 * transform( HT_ang_U2D1, j );
/// frame4 is computed from frame3 by walking along frame3's z-axis by dis_U1D1_samples( k );
/// frame4.point() is the coordinate for atom 4.
/// frame5 is the product: frame4 * transform( HT_tor_U2D2, l );
/// frame6 is the product: frame5 * transform( HT_ang_U1D2, m );
/// frame6.point() is the coordinate for atom 5.
/// NOTE: transform( HT_ang_U1D2, m ) is the product of the z-axis
/// rotation by ang_U1D2_sample_[ m ] and the transformation HT along the
/// z-axis by the distance between atoms D1 and D2.  The ExternalGeomSampler
/// pre-multiplies these matricies.
/// frame7 is the product: frame6 * transform( HT_tor_U1D3, n );
/// frame7.point() is the coordinate for atom 6.
/// NOTE: transform( HT_tor_U1D3, n ) is the product of three HTs:
/// the HT representing the z-axis rotation by tor_U1D3_sample_[ n ],
/// the HT representing the x-axis rotation by ( -1 * ( 180 - ang_D1D2D3 ) )
/// and the HT representing the z-axis translation by the
/// distance between atoms D2 and D3.  The ExternalGeomSampler
/// pre-multiplies these matrices.

class ExternalGeomSampler : public utility::pointer::ReferenceCount {
public:
	typedef utility::pointer::ReferenceCount parent;
	typedef core::Size Size;
	typedef core::Real Real;
	typedef numeric::HomogeneousTransform< Real > HTReal;

public:
	virtual ~ExternalGeomSampler();

	ExternalGeomSampler();
	ExternalGeomSampler( ExternalGeomSampler const & );

	ExternalGeomSampler &
	operator = ( ExternalGeomSampler const & rhs );


public:
	/// Intialization routines
	/// all initialization routines must be invoked before building may begin.
	/// Once each has been called, the function precompute_transforms() must
	/// be called.

	void set_tor_U3D1_samples( utility::vector1< Real > const & tor_U3D1_samples );
	void set_dis_U1D1_samples( utility::vector1< Real > const & dis_U1D1_samples );
	void set_ang_U2D1_samples( utility::vector1< Real > const & ang_U2D1_samples );
	void set_ang_U1D2_samples( utility::vector1< Real > const & ang_U1D2_samples );
	void set_tor_U2D2_samples( utility::vector1< Real > const & tor_U2D2_samples );
	void set_tor_U1D3_samples( utility::vector1< Real > const & tor_U1D3_samples );

	void set_tor_U3D1_samples( std::list< Real > const & tor_U3D1_samples );
	void set_dis_U1D1_samples( std::list< Real > const & dis_U1D1_samples );
	void set_ang_U2D1_samples( std::list< Real > const & ang_U2D1_samples );
	void set_ang_U1D2_samples( std::list< Real > const & ang_U1D2_samples );
	void set_tor_U2D2_samples( std::list< Real > const & tor_U2D2_samples );
	void set_tor_U1D3_samples( std::list< Real > const & tor_U1D3_samples );


	void set_dis_D1D2(    Real distance );
	void set_dis_D2D3(    Real distance );
	void set_ang_D1D2D3( Real ang_in_degrees );

	/// @brief Must be called after the samples are set, and the
	/// internal geometry of the three downstream coordinates
	/// (point 4, 5, and 6) are described.  Does nothing if the transforms
	/// are up to date.
	void
	precompute_transforms();

	/// @brief convenience function if one wants to change the meaning of
	/// upstream and downstream
	void
	flip_upstream_downstream_samples();

public:
	/// Accessors

	Size n_tor_U3D1_samples() const { return tor_U3D1_samples_.size(); }
	Size n_dis_U1D1_samples() const { return dis_U1D1_samples_.size(); }
	Size n_ang_U2D1_samples() const { return ang_U2D1_samples_.size(); }
	Size n_ang_U1D2_samples() const { return ang_U1D2_samples_.size(); }
	Size n_tor_U2D2_samples() const { return tor_U2D2_samples_.size(); }
	Size n_tor_U1D3_samples() const { return tor_U1D3_samples_.size(); }

	utility::vector1< Real > const & tor_U3D1_samples() const { return tor_U3D1_samples_; }
	utility::vector1< Real > const & dis_U1D1_samples() const { return dis_U1D1_samples_; }
	utility::vector1< Real > const & ang_U2D1_samples() const { return ang_U2D1_samples_; }
	utility::vector1< Real > const & ang_U1D2_samples() const { return ang_U1D2_samples_; }
	utility::vector1< Real > const & tor_U2D2_samples() const { return tor_U2D2_samples_; }
	utility::vector1< Real > const & tor_U1D3_samples() const { return tor_U1D3_samples_; }

public:
	bool transforms_uptodate() const { return transforms_uptodate_; }


	HTReal const &
	transform( ExternalTransform id, Size which_state ) const {
		assert( transforms_uptodate_ );
		return transforms_[ id ][ which_state ];
	}

private:
	utility::vector1< Real > dis_U1D1_samples_;
	utility::vector1< Real > ang_U2D1_samples_;
	utility::vector1< Real > tor_U3D1_samples_;
	utility::vector1< Real > ang_U1D2_samples_;
	utility::vector1< Real > tor_U2D2_samples_;
	utility::vector1< Real > tor_U1D3_samples_;

	Real dis_D1D2_;
	Real dis_D2D3_;
	Real ang_D1D2D3_;

	bool transforms_uptodate_;

	utility::vector1< utility::vector1< HTReal > > transforms_;

};

}
}
}

#endif
