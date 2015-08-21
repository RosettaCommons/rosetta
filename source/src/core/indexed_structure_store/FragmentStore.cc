// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file
/// @brief
/// @author Alex Ford <fordas@uw.edu>

#include <boost/algorithm/string/join.hpp>
#include <boost/format.hpp>
#include <core/indexed_structure_store/FragmentStore.hh>
#include <numeric/xyzVector.hh>

namespace core
{
namespace indexed_structure_store
{
std::ostream& operator<<(std::ostream& os, const FragmentSpecification& s)
{
	boost::format f("FragmentSpecification(fragment_length=%s, fragment_atoms=[%s])");
	os << f % s.fragment_length % boost::algorithm::join(s.fragment_atoms, ", ");
	return os;
}

FragmentStore::FragmentStore(
	FragmentSpecification fragment_specification,
	numeric::Size num_fragments) :
	fragment_specification(fragment_specification),
	fragment_threshold_distances(num_fragments),
	fragment_coordinates(num_fragments * fragment_specification.coordinates_per_fragment() )
{
}

void FragmentStore::resize(numeric::Size num_fragments)
{
	fragment_threshold_distances.resize(num_fragments);
	fragment_coordinates.resize(num_fragments * fragment_specification.coordinates_per_fragment());
}

void FragmentStore::add_threshold_distance_allFrag(numeric::Real distance)
{
	for ( numeric::Size ii=0; ii< fragment_threshold_distances.size(); ++ii ) {
		fragment_threshold_distances[ii]=distance;
	}
}

}
}
