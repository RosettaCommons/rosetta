// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/geometry/BoundingBox.hh
/// @brief  3d axis aligned bounding box
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_numeric_geometry_BoundingBox_hh
#define INCLUDED_numeric_geometry_BoundingBox_hh


// unit headers
#include <numeric/geometry/BoundingBox.fwd.hh>


namespace numeric {
namespace geometry {


/// @brief 3d axis aligned bounding box class
/// @details template type must be 3d and have .x(), .y(), .z() accessors
template< typename T >
class BoundingBox {

public: // types

	/// @brief 3d point position
	typedef T PointPosition;


public: // construct/destruct

	/// @brief default constructor
	/// @warning no initialization of corners for speed, make sure you reset()
	/// @warning or otherwise set corners before adding points
	inline
	BoundingBox() {}

	/// @brief point constructor
	inline
	BoundingBox(
		PointPosition const & pp
	) : lower_( pp ),
		upper_( pp )
	{}

	/// @brief corner constructor
	inline
	BoundingBox(
		PointPosition const & lower,
		PointPosition const & upper
	) : lower_( lower ),
		upper_( upper )
	{}

	/// @brief copy constructor
	inline
	BoundingBox(
		BoundingBox const & bb
	) : lower_( bb.lower_ ),
		upper_( bb.upper_ )
	{}

	/// @brief default destructor
	inline
	~BoundingBox() {}


public: // assignment

	/// @brief copy assignment
	inline
	BoundingBox &
	operator =(
		BoundingBox const & bb
	)
	{
		if ( this != &bb ) {
			lower_ = bb.lower_;
			upper_ = bb.upper_;
		}
		return *this;
	}


public: // box management

	/// @brief add a point, expands bounds if necessary
	inline
	void
	add(
		PointPosition const & pp
	)
	{
		lower_.min( pp );
		upper_.max( pp );
	}

	/// @brief get lower corner
	inline
	PointPosition const &
	lower() const
	{
		return lower_;
	}

	/// @brief get upper corner
	inline
	PointPosition const &
	upper() const
	{
		return upper_;
	}

	/// @brief set lower corner
	inline
	void
	set_lower(
		PointPosition const & p
	)
	{
		lower_ = p;
	}

	/// @brief set upper corner
	inline
	void
	set_upper(
		PointPosition const & p
	)
	{
		upper_ = p;
	}

	/// @brief reset corners
	inline
	void
	reset(
		PointPosition const & p = PointPosition()
	)
	{
		lower_ = p;
		upper_ = p;
	}

	/// @brief expand box corners (additive)
	template< typename U >
	inline
	void
	expand(
		U const & scalar
	)
	{
		lower_ -= scalar;
		upper_ += scalar;
	}

	// @brief contract box corners (subtractive)
	template< typename U >
	inline
	void
	contract(
		U const & scalar
	)
	{
		lower_ += scalar;
		upper_ -= scalar;
	}

	/// @brief translate bounding box
	inline
	void
	translate(
		PointPosition const & t
	)
	{
		lower_ += t;
		upper_ += t;
	}


public: // box query

	/// @brief intersects another bounding box?
	inline
	bool
	intersects(
		BoundingBox const & bb
	) const
	{
		return !( lower_.x() > bb.upper_.x() || bb.lower_.x() > upper_.x() ||
			lower_.y() > bb.upper_.y() || bb.lower_.y() > upper_.y() ||
			lower_.z() > bb.upper_.z() || bb.lower_.z() > upper_.z() );
	}

	/// @brief is point contained within this bounding box?
	template< typename U >
	inline
	bool
	contains(
		U const & x,
		U const & y,
		U const & z
	) const
	{
		return lower_.x() <= x && lower_.y() <= y && lower_.z() <= z &&
			x <= upper_.x() && y <= upper_.y() && z <= upper_.z();
	}

	/// @brief is point contained within this bounding box?
	inline
	bool
	contains(
		PointPosition const & p
	) const
	{
		return contains( p.x(), p.y(), p.z() );
	}


private: // data

	PointPosition lower_; // lower corner
	PointPosition upper_; // upper corner

};


}
}


#endif /*INCLUDED_numeric_geometry_BoundingBox_HH*/
