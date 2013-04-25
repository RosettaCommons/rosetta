// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file WrappedPrimitive
/// @brief silly little class to wrap primitive types in an OP.
/// @author James Thompson

#ifndef INCLUDED_numeric_kdtree_WrappedType_hh
#define INCLUDED_numeric_kdtree_WrappedType_hh

#include <numeric/types.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>


namespace numeric {
namespace kdtree {

template< typename T >
class WrappedPrimitive : public utility::pointer::ReferenceCount {
public:
	WrappedPrimitive( T val ) : val_( val ) {}

	T val() const {
		return val_;
	}

	void val( T new_val ) {
		new_val = val_;
	}

private:
	T val_;
};

} // kdtree
} // numeric

#endif
