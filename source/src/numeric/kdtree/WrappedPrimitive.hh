// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/kdtree/WrappedPrimitive
/// @brief silly little class to wrap primitive types in an OP.
/// @author James Thompson

#ifndef INCLUDED_numeric_kdtree_WrappedPrimitive_hh
#define INCLUDED_numeric_kdtree_WrappedPrimitive_hh

#include <utility/pointer/ReferenceCount.hh>


namespace numeric {
namespace kdtree {

// Requirements: Class T must have a default constructor and an assignment operator.
template< typename T >
class WrappedPrimitive : public utility::pointer::ReferenceCount {
public:
	WrappedPrimitive() : val_() {}
	WrappedPrimitive( T const & val ) : val_( val ) {}

	T & val() {
		return val_;
	}

	T const & val() const {
		return val_;
	}

	void val( T const & new_val ) {
		val_ = new_val;
	}

private:
	T val_;
};

} // kdtree
} // numeric

#endif
