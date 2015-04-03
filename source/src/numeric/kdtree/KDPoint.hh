// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file numeric/kdtree/KDPoint.hh
/// @brief utility functions for kd-tree. See kdtree.hh for more information.
/// @author James Thompson
//
#ifndef INCLUDED_numeric_kdtree_KDPoint_hh
#define INCLUDED_numeric_kdtree_KDPoint_hh

#include <numeric/types.hh>
#include <numeric/kdtree/KDPoint.fwd.hh>

#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace numeric {
namespace kdtree {

class KDPoint : public utility::pointer::ReferenceCount {

public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~KDPoint();
	KDPoint();

	KDPoint(
		utility::vector1< numeric::Real > location
	);

	KDPoint(
		utility::vector1< numeric::Real > location,
		utility::pointer::ReferenceCountOP data
	);

	KDPoint(
		utility::vector1< numeric::Real > location,
		numeric::Real distance
	);

	KDPoint(
		utility::vector1< numeric::Real > location,
		utility::pointer::ReferenceCountOP data,
		numeric::Real distance
	);

	/// @brief copy constructor
	KDPoint( KDPoint const & src );

	/// @brief Returns a const reference to the location of this point in
	/// k-space.
	utility::vector1< numeric::Real > const & const_location() const;

	/// @brief Returns the number of dimensions for the space in which this point
	/// lives.
	numeric::Size size() const;

	/// @brief Returns to the location of this point in k-space.
	utility::vector1< numeric::Real > location() const;

	/// @brief Returns the data associated with this class.
	utility::pointer::ReferenceCountOP data() const;

	/// @brief getter for distance() from this point to an arbitrary point
	/// in kd-space.
	numeric::Real distance() const;

	/// @brief Sets the distance to an arbitrary point.
	void distance( numeric::Real dist );

	/// @brief sets the location of this point in k-space.
	void location( utility::vector1< numeric::Real > dat );

	/// @brief Assignment operator for KDPoint class.
	KDPoint & operator = ( KDPoint const & src );

	/// @brief Equality operator. Compares location() and distance(),
	/// ignores data().
	bool operator == ( KDPoint const & other ) const;

	/// @brief Prints the definition of this object to the given ostream.
	void show( std::ostream & out ) const;

	/// @brief Returns a stringified version of this object.
	std::string to_string() const;

	/// @brief Reads the definition for this object from the given istream.
	void read_data( std::istream & in );

private:
	// location in k-space of this point
	utility::vector1< numeric::Real > location_;

	// data associated with this point
	utility::pointer::ReferenceCountOP data_;

	// distance to an arbitrary point, intended for use in nearest-neighbor
	// search.
	numeric::Real distance_;
};

// simple class for comparing points in k-space by the value in the k-th
// dimension.
class CompareKDPoints {
public:
	CompareKDPoints() : idx_(1) {}
	CompareKDPoints( numeric::Size new_idx ) : idx_(new_idx) {}

	bool operator() (
		KDPointOP x, KDPointOP y
	) {
		utility::vector1< numeric::Real > const & x_dat( x->const_location() );
		utility::vector1< numeric::Real > const & y_dat( y->const_location() );
		assert( x->location().size() == y->location().size() );
		assert( idx_ <= x->location().size() );
		return ( x_dat[idx_] < y_dat[idx_] );
	}

private:
	numeric::Size idx_;
};

} // kdtree
} // numeric

#endif
