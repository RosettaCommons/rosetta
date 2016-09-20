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

#include <boost/algorithm/string/join.hpp>
#include <boost/format.hpp>
#include <core/indexed_structure_store/FragmentStore.hh>
#include <core/indexed_structure_store/FragmentLookup.hh>

namespace core
{
namespace indexed_structure_store
{
std::ostream& operator<<(std::ostream& os, const FragmentSpecification& s)
{
	os << "FragmentSpecification(fragment_length=" << s.fragment_length
		<< ", fragment_atoms=[" << boost::algorithm::join(s.fragment_atoms, ", ") << "])";
	return os;
}

FragmentStore::FragmentStore(
	FragmentSpecification fragment_specification,
	numeric::Size num_fragments) :
	fragment_specification(fragment_specification),
	fragment_threshold_distances(num_fragments),
	fragment_coordinates(num_fragments * fragment_specification.coordinates_per_fragment() )
{
	num_fragments_ = num_fragments;
	fragLookupOP_=nullptr;
}

void FragmentStore::resize(numeric::Size num_fragments)
{
	num_fragments_ = num_fragments;
	fragment_threshold_distances.resize(num_fragments);
	fragment_coordinates.resize(num_fragments * fragment_specification.coordinates_per_fragment());
}

void FragmentStore::add_threshold_distance_allFrag(numeric::Real distance)
{
	for ( double & fragment_threshold_distance : fragment_threshold_distances ) {
		fragment_threshold_distance=distance;
	}
}

void FragmentStore::generate_subset_fragment_store(std::vector<numeric::Size> residues, std::string name){
	numeric::Size tmp_fragment_length = residues.size();
	FragmentSpecification tmp_fragment_spec = FragmentSpecification(tmp_fragment_length,fragment_specification.fragment_atoms);
	FragmentStoreOP tmp_fragment_store = FragmentStoreOP(new FragmentStore(tmp_fragment_spec,num_fragments_));
	for ( numeric::Size ii=0; ii<num_fragments_; ++ii ) {//num fragments
		Size tmp_frag_index=0;
		for ( numeric::Size jj=0; jj<fragment_specification.coordinates_per_fragment(); ++jj ) {
			if ( std::find(residues.begin(),residues.end(),jj)!= residues.end() ) { //only copy the residue that are needed
				for ( numeric::Size kk=0; kk<fragment_specification.fragment_atoms.size(); ++kk ) {//currently always 1
					for ( numeric::Size ll=0; ll<3; ++ll ) {//number of residues per atom
						Size start_index1 = ii*fragment_specification.coordinates_per_fragment()+jj*fragment_specification.fragment_atoms.size()+kk;
						Size start_index2 = ll;
						Size end_index1 = ii*tmp_fragment_spec.coordinates_per_fragment()+tmp_frag_index*fragment_specification.fragment_atoms.size()+kk;
						Size end_index2 = ll;
						numeric::Real start_position=fragment_coordinates[start_index1][start_index2];
						tmp_fragment_store->fragment_coordinates[end_index1][end_index2]=start_position;
						//numeric::Real *start_position=&fragment_coordinates[start_index1][start_index2];
						//numeric::Real *end_position=&tmp_fragment_store->fragment_coordinates[end_index1][end_index2];
						//end_position=start_position;
					}
				}
				tmp_frag_index++;
			}
		}
	}
	tmp_fragment_store->add_threshold_distance_allFrag(0);
	fragmentStore_groups[name]=tmp_fragment_store;
}


//done so that the fragmentLookup is cached when calculated.
FragmentLookupOP FragmentStore::get_fragmentLookup(){
	if ( fragLookupOP_==nullptr ) {
		fragLookupOP_ = FragmentLookupOP(new FragmentLookup(this->get_self_ptr()));
	}
	return(fragLookupOP_);
}

//get fragment coordinates
std::vector<numeric::xyzVector<numeric::Real> > FragmentStore::get_fragment_coordinates(numeric::Size position){
	std::vector<numeric::xyzVector<numeric::Real> > tmp_coordinates;
	Size start_position = position*fragment_specification.coordinates_per_fragment()*fragment_specification.fragment_atoms.size();
	for ( numeric::Size ii=0; ii<fragment_specification.coordinates_per_fragment(); ++ii ) {
		for ( numeric::Size jj=0; jj<fragment_specification.fragment_atoms.size(); ++jj ) {
			Size db_coord = start_position+ii*fragment_specification.fragment_atoms.size()+jj;
			//Size frag_coord = ii*fragment_specification.fragment_atoms.size()+jj;
			tmp_coordinates.push_back(fragment_coordinates[db_coord]);
		}
	}
	return(tmp_coordinates);
}


}
}
