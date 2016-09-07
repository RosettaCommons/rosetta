// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/kdtree/KDPointList.hh
/// @brief definition for a class that keeps track of the best N KDPoint
/// objects by distance.
/// @author James Thompson

#ifndef INCLUDED_numeric_kdtree_KDPointList_hh
#define INCLUDED_numeric_kdtree_KDPointList_hh

#include <numeric/types.hh>
#include <numeric/kdtree/KDPoint.hh>

#include <utility/vector1.hh>

namespace numeric {
namespace kdtree {

struct KDPoint_MinDist {
	bool operator()( KDPointOP const & a, KDPointOP const & b );
};

/// @brief Class for keeping track of the closest N KDPoint
/// objects by distance.
class KDPointList : public utility::pointer::ReferenceCount {
public: // iterators
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~KDPointList() override;

	typedef utility::vector1< KDPointOP >::iterator iterator;
	typedef utility::vector1< KDPointOP >::const_iterator const_iterator;
	const_iterator begin() const;
	const_iterator end() const;
	iterator begin();
	iterator end();

public:
	KDPointList( numeric::Size const n_to_keep );

	void insert( KDPointOP pt );

	numeric::Real worst_distance() const;

	KDPointOP operator[]( numeric::Size const pos ) const;

	numeric::Size size() const;

	numeric::Size max_values() const;

	numeric::Real distance_cutoff() const;
	void distance_cutoff( numeric::Real const cutoff );

	utility::vector1< KDPointOP > sorted_values();

	/// @brief merge another KDPointList with this KDPointList.
	/// This calls insert which is a little slow, and is a candidate for
	/// optimization if the insert() method shows up in profiling.
	void merge( KDPointList const & other );

	void show( std::ostream & out ) const;

private:
	void update_heap_();

	void update_size_();

private:
	numeric::Size max_vals_;
	numeric::Real distance_cutoff_;
	utility::vector1< KDPointOP > container_;
};

} // kdtree
} // numeric

#endif
