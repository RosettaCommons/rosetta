// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/numeric/kdtree/KDPoint.cc
/// @brief
/// @author James Thompson


#include <numeric/types.hh>
#include <numeric/kdtree/WrappedReal.hh>
#include <numeric/kdtree/KDPoint.hh>

#include <utility/vector1.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>

#include <ObjexxFCL/string.functions.hh>

namespace numeric {
namespace kdtree {

/// @details Auto-generated virtual destructor
KDPoint::~KDPoint() {}

KDPoint::KDPoint() {}

KDPoint::KDPoint(
	utility::vector1< numeric::Real > location
) :
	location_( location ), data_( /* NULL */ ), distance_( 0.0 )
{}

KDPoint::KDPoint(
	utility::vector1< numeric::Real > location,
	utility::pointer::ReferenceCountOP data
) :
	location_( location ), data_( data ), distance_( 0.0 )
{}

KDPoint::KDPoint(
	utility::vector1< numeric::Real > location,
	numeric::Real distance
) :
	location_( location ), data_( /* NULL */ ), distance_( distance )
{}

KDPoint::KDPoint(
	utility::vector1< numeric::Real > location,
	utility::pointer::ReferenceCountOP data,
	numeric::Real distance
) :
	location_( location ), data_( data ), distance_( distance )
{}

KDPoint::KDPoint( KDPoint const & src ) : ReferenceCount() {
	// use assignment operator
	*this = src;
}

/// @brief Returns a const reference to the location of this point in
///	k-space.
utility::vector1< numeric::Real > const & KDPoint::const_location() const {
	return location_;
}

/// @brief Returns the number of dimensions for the space in which this point
/// lives.
numeric::Size KDPoint::size() const {
	return location_.size();
}

/// @brief Returns to the location of this point in k-space.
utility::vector1< numeric::Real > KDPoint::location() const {
	return location_;
}

utility::pointer::ReferenceCountOP KDPoint::data() const {
	return data_;
}

numeric::Real KDPoint::distance() const {
	return distance_;
}

/// @brief Sets the distance to an arbitrary point.
void KDPoint::distance( numeric::Real dist ) {
	distance_ = dist;
}

/// @brief sets the location of this point in k-space.
void KDPoint::location( utility::vector1< numeric::Real > dat ) {
	location_ = dat;
}

/// @brief Assignment operator for KDPoint class.
KDPoint & KDPoint::operator = ( KDPoint const & src ) {
	location_   = src.location();
	distance_   = src.distance();
	if ( src.data() ) {
		data_       = src.data();
	}
	return *this;
}

/// @brief Equality operator. Compares location() and distance(),
/// ignores data().
bool KDPoint::operator == ( KDPoint const & other ) const {
	return (
		location()   == other.location() &&
		distance() == other.distance()
	);
}

void KDPoint::show( std::ostream & out ) const {
	out << to_string() << " distance = " << distance();
}

/// should look like:
/// KDPOINT ndim val1 val2 val3 DATA <data definition>
std::string KDPoint::to_string() const {
	using ObjexxFCL::string_of;
	std::string retval("KDPOINT ");

	retval += string_of( size() ) + " ";
	for ( utility::vector1< Real >::const_iterator
				it = location_.begin(), end = location_.end() - 1;
				it != end; ++it
	) {
		retval += string_of( *it );
		retval += " ";
	}

	// do something smarter with data here!
	if ( data() ) {
		retval += "DATA ";
		WrappedRealOP val = utility::pointer::dynamic_pointer_cast< WrappedReal > ( data() );
		retval += string_of( val->val() );
	}

	return retval;
}

void KDPoint::read_data( std::istream & in ) {
	std::string token;
	in >> token;
	assert( token == "KDPOINT" );
	Size ndim;
	in >> ndim;

	location_.resize( ndim, 0.0 );
	for ( Size ii = 1; ii <= ndim; ++ii ) {
		Real val;
		in >> val;
		location_[ ii ] = val;
	}
	// do something smarter with data here!
	in >> token;
	assert( token == "DATA" );
	if ( in.good() ) {
		Real val;
		in >> val;
		WrappedRealOP dat( new WrappedReal( val ) );
		data_ = dat;
	}
}

} // kdtree
} // numeric
