// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/kdtree/HyperRectangle.hh
/// @brief
/// @author James Thompson
//
#ifndef INCLUDED_numeric_kdtree_HyperRectangle_hh
#define INCLUDED_numeric_kdtree_HyperRectangle_hh

#include <numeric/types.hh>
#include <numeric/kdtree/HyperRectangle.fwd.hh>

#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace numeric {
namespace kdtree {

class HyperRectangle : public utility::pointer::ReferenceCount {

public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~HyperRectangle();
	HyperRectangle();

	HyperRectangle(
		utility::vector1< numeric::Real > upper,
		utility::vector1< numeric::Real > lower
	);

	HyperRectangle(
		utility::vector1< utility::vector1< Real > > const & pts
	);

	utility::vector1< numeric::Real > upper() const;
	utility::vector1< numeric::Real > lower() const;

	numeric::Size ndim() const;

	HyperRectangle & operator = ( HyperRectangle const & src );

	void extend( utility::vector1< numeric::Real > const & pt );

	void show( std::ostream & out ) const;

private:
	utility::vector1< numeric::Real > upper_, lower_;
};

} // kdtree
} // numeric

#endif
