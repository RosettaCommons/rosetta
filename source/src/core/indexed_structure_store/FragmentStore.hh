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

#ifndef INCLUDED_core_indexed_structure_store_FragmentStore_hh
#define INCLUDED_core_indexed_structure_store_FragmentStore_hh

#include <utility/pointer/ReferenceCount.hh>

#include <numeric/types.hh>

#include <core/indexed_structure_store/FragmentStore.fwd.hh>

#include <vector>
#include <string>
#include <iostream>

namespace core
{
namespace indexed_structure_store
{

struct FragmentSpecification
{
	FragmentSpecification() :
		fragment_length(0),
		fragment_atoms(0)
	{
		
	}

	FragmentSpecification(numeric::Size fragment_length, std::vector<std::string> fragment_atoms) :
		fragment_length(fragment_length),
		fragment_atoms(fragment_atoms)
	{
		
	}

	numeric::Size fragment_length;
	std::vector<std::string> fragment_atoms;

	numeric::Size coordinates_per_fragment() const { return fragment_atoms.size() * fragment_length; }

	friend std::ostream& operator<<(std::ostream& os, const FragmentSpecification& s);
};

class FragmentStore : public utility::pointer::ReferenceCount
{
public:
	// @brief Basic structure store, holds a collection of structure and associated residue entries.
	FragmentStore(FragmentSpecification fragment_specification, numeric::Size num_fragments = 0);

	void resize(numeric::Size num_fragments);

  FragmentSpecification fragment_specification;
	std::vector<numeric::Real> fragment_threshold_distances;
	std::vector< numeric::xyzVector<numeric::Real> > fragment_coordinates;
};

}
}
#endif
