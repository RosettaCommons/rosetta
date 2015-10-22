// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   /numeric/geometry/hashing/SixDHasher.fwd.hh
/// @brief  Declaration for classes in 6D hasher
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini
/// @author Ken Jung (kenjung@uw.edu), adding radial lookup


#ifndef INCLUDED_numeric_geometry_hashing_SixDHasher_hh
#define INCLUDED_numeric_geometry_hashing_SixDHasher_hh

// Unit headers
#include <numeric/geometry/hashing/SixDHasher.fwd.hh>


// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/geometry/BoundingBox.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

/// Boost headers
#include <boost/unordered_map.hpp>

/// C++ headers
#include <vector>

//#include <protocols/match/Hit.fwd.hh>
#include <utility/fixedsizearray1.hh>


namespace numeric {
namespace geometry {
namespace hashing {

/// @brief Small hashing struct with no private data; simply an algorithm to turn a 64-bit
/// representation of a 6d voxel into an integer.
struct bin_index_hasher
{
	/// @brief functor used by boost (and sgi's stl) hash classes.
	boost::uint64_t
	operator() ( boost::uint64_t bin_index ) const {
		boost::uint64_t const P = 999987;   // largest prime  number less than max int value coded by std::size_t int (32 bits)
		return  bin_index % P;
	}

};

/// @brief Returns a list of offsets corresponding to the bins in a hypershell with radius x
class SixDOffsetTree {
public:
	typedef numeric::Real                               Real;
	typedef numeric::Size                               Size;
	typedef platform::SSize                              SSize;

	SixDOffsetTree();

	// returns only offsets within bounds
	std::vector< SBin6D > lookup( Size radius, const Bin6D & center, const Bin6D & bounds ) const;

	//don't put this in the constructor because it will add overhead to anyone using the 6dhasher without the radial tree
	void init( Size max_radius );
	Size sum_radius( SBin6D & input, Size range = 6);


private:
	bool insert( SBin6D & input, Size depth = 1, Size caller = 0);

	// Holds the offsets
	// instead of pointers we use integer offsets
	// format: data_[radius][random idx][depth]
	// either points to the next depth or 0/1 if depth=6
	std::vector < std:: vector < boost::unordered_map < SSize, Size > > > data_;
};


/// @brief Bin the six degrees of freedom that describe the downstream geometry of a hit.
/// These degrees of freedom are, in order, the x, y and z coordinates of orientation atom3,
/// and the phi, psi, and theta euler angles that describe the orientation of the coordinate
/// frame at orientation atom 3.  The binner is responsible for maintaining the lower corner
/// of the 6-d space -- the first two Euler angles  wrap at 360; the third Euler
/// angle, theta, does not wrap in the same way.  See the comments for the bin6 method.
class SixDCoordinateBinner : public utility::pointer::ReferenceCount
{
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~SixDCoordinateBinner();
	typedef numeric::Real                               Real;
	typedef numeric::Size                               Size;
	typedef numeric::xyzVector< numeric::Real >         Vector;
	typedef numeric::geometry::BoundingBox< Vector > BoundingBox;

private:
	SixDCoordinateBinner();

public:
	SixDCoordinateBinner(
		BoundingBox const & bounding_box,
		Size3 const & euler_offsets,
		utility::fixedsizearray1< Real, 6 > bin_widths
	);

	inline
	bool
	contains( Real6 const & point ) const {
		return bounding_box_.contains( Vector( point[1], point[2], point[3] ));
	}

	void tree_init( Size max_radius ) {
		offset_tree_.init( max_radius );
	}


	/// @brief Construct the discrete representation of a six-dimensional vector
	/// of reals representing an xyz coordinate in its first three dimensions and
	/// a set of three Euler angles in the last three dimensions.
	/// Precondition: The xyz coordinate must be inside the bounding box of this binner.
	/// Precondition: The euler angles should be in degrees; the first two should
	/// be in the range between 0 and 360, the third should be in the range from
	/// 0 to 180.
	Bin6D
	bin6( Real6 const & values ) const;

	/// @brief Determine halfbin index for a point in 6D
	/// i.e., how far from the lower corner of the point's containing 6D voxel
	/// is the point -- is it more than halfway to the next 6D voxel?
	/// Each dimension returned will hold a 0 or a 1.
	Bin6D
	halfbin6( Real6 const & values ) const;

	/// @brief functor used by boost (and sgi's stl) hash classes.
	boost::uint64_t
	bin_index( Bin6D const & bin ) const {

		assert( bin[ 1 ] < dimsizes_[ 1 ] );
		assert( bin[ 2 ] < dimsizes_[ 2 ] );
		assert( bin[ 3 ] < dimsizes_[ 3 ] );
		assert( bin[ 4 ] < dimsizes_[ 4 ] );
		assert( bin[ 5 ] < dimsizes_[ 5 ] );
		assert( bin[ 6 ] < dimsizes_[ 6 ] );

		boost::uint64_t const A =
			bin[ 1 ] * dimprods_[ 1 ] +
			bin[ 2 ] * dimprods_[ 2 ] +
			bin[ 3 ] * dimprods_[ 3 ] +
			bin[ 4 ] * dimprods_[ 4 ] +
			bin[ 5 ] * dimprods_[ 5 ] +
			bin[ 6 ] * dimprods_[ 6 ];

		return A;

	}

	/// @brief compute the bin index (64-bit int) for a 6D point.
	boost::uint64_t
	bin_index( Real6 const & values ) const {
		Bin6D bin = bin6( values );
		return bin_index( bin );
	}

	// @brief compute the list of bin indices of a hypershell
	// Bounding box overflow pruned here
	std::vector < boost::uint64_t >
	radial_bin_index( numeric::Size radius, Real6 const & center ) const;

	Bin6D
	bin_from_index( boost::uint64_t index ) const {
		Bin6D bin;
		for ( Size ii = 1; ii <= 6; ++ii ) {
			bin[ ii ] = index / dimprods_[ ii ];
			index = index % dimprods_[ ii ];
			//std::cout << "bin:  " << bin[ii] << std::endl;
		}
		return bin;
	}

	Real6
	bin_center_point( Bin6D const & bin ) const;

	BoundingBox
	bounding_volume_from_index( boost::uint64_t index ) const {
		Bin6D bin = bin_from_index( index );
		Vector lower(
			bin[ 1 ] * bin_widths_[ 1 ] + bounding_box_.lower()( 1 ),
			bin[ 2 ] * bin_widths_[ 2 ] + bounding_box_.lower()( 2 ),
			bin[ 3 ] * bin_widths_[ 3 ] + bounding_box_.lower()( 3 ) );
		Vector upper( bin_widths_[ 1 ], bin_widths_[ 2 ], bin_widths_[ 3 ] );
		upper += lower;
		return BoundingBox( lower, upper );

	}

	BoundingBox const & bounding_box() const { return bounding_box_; }
	Real3 const & euler_offsets() const { return euler_offsets_; }
	Real6 const & bin_widths() const { return bin_widths_; }
	Real6 const & halfbin_widths() const { return halfbin_widths_; }

	Size6 const &
	dimsizes() const {
		return dimsizes_;
	}

private:
	/// When using a particular set of offsets, make sure that the euler angles wrap so that
	/// neighboring points in 6D actually fall into the same bins.  The logic here is tested
	/// thuroughly in Hasher.cxxtest.hh
	Real3 wrap_euler_angles( Real6 const & values ) const;

private:
	BoundingBox bounding_box_;
	Size6 dimsizes_;
	Size6 dimprods_;
	Real6 bin_widths_;
	Real6 halfbin_widths_;
	Real3 euler_offsets_;
	SixDOffsetTree offset_tree_;
};


}//hashing
}//geometry
}//numeric

#endif
