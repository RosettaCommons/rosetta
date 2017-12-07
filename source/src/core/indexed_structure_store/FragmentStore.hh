// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file
/// @brief
/// @author Alex Ford <fordas@uw.edu>

#ifndef INCLUDED_core_indexed_structure_store_FragmentStore_hh
#define INCLUDED_core_indexed_structure_store_FragmentStore_hh

#include <utility>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>

#include <numeric/xyzVector.fwd.hh>
#include <numeric/types.hh>

#include <core/indexed_structure_store/FragmentStore.fwd.hh>
#include <core/indexed_structure_store/FragmentLookup.fwd.hh>

#include <vector>
#include <map>
#include <set>
#include <iosfwd>

#include <ctime>

namespace core {
namespace indexed_structure_store {

struct FragmentSpecification
{
	FragmentSpecification() :
		fragment_length(0),
		fragment_atoms(0)
	{

	}

	FragmentSpecification(numeric::Size fragment_length, std::vector<std::string> fragment_atoms) :
		fragment_length(fragment_length),
		fragment_atoms(std::move(fragment_atoms))
	{

	}

	numeric::Size fragment_length;
	std::vector<std::string> fragment_atoms;

	numeric::Size coordinates_per_fragment() const { return fragment_atoms.size() * fragment_length; }

	friend std::ostream& operator<<(std::ostream& os, const FragmentSpecification& s);
};

class FragmentStore : public utility::pointer::ReferenceCount , public utility::pointer::enable_shared_from_this< FragmentStore >

{
public:
	// @brief Basic structure store, holds a collection of structure and associated residue entries.
	FragmentStore(FragmentSpecification fragment_specification, numeric::Size num_fragments = 0);

	/// self pointer
	inline FragmentStoreOP get_self_ptr() { return shared_from_this(); }


	void resize(numeric::Size num_fragments);
	void add_threshold_distance_allFrag(numeric::Real distance);

	template <class T> void inline_vector_delete(std::vector<T>& data_vector, const std::vector<bool>& delete_vector){
		typename std::vector<T>::iterator read_position = data_vector.begin();
		typename std::vector<T>::iterator write_position = data_vector.begin();
		typename std::vector<T>::iterator end_position = data_vector.end();
		//for (std::_Bit_const_iterator::const_reference itr : delete_vector) {
		for ( bool const itr : delete_vector ) {
			if ( itr==true ) {
				read_position++;
			} else {
				*write_position = *read_position;
				read_position++;
				write_position++;
			}
		}
		data_vector.erase(write_position,end_position);
	}

	std::set<std::string> get_homologs();
	void delete_homologs();

	void generate_residue_subset_fragment_store(std::vector<numeric::Size> residues);

	FragmentLookupOP get_fragmentLookup();

	std::vector<numeric::xyzVector<numeric::Real> > get_fragment_coordinates(numeric::Size position);


	FragmentSpecification fragment_specification;
	numeric::Size num_fragments_;
	numeric::Size hash_id;
	std::vector<numeric::Real> fragment_threshold_distances;
	std::vector<numeric::xyzVector<numeric::Real> > fragment_coordinates;
	std::map<std::string, std::vector<numeric::Size> > int64_groups;
	std::map<std::string, std::vector<numeric::Real> > real_groups;
	std::map<std::string, std::vector<std::vector<numeric::Real> > >realVector_groups;
	std::map<std::string, std::vector<std::string> > string_groups;
	FragmentStoreOP residue_subset_fragment_store;

private:
	FragmentLookupOP fragLookupOP_; //stored here to cache the removal of center of mass
};

}
}
#endif
