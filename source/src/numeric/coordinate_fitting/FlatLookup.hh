// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <vector>
#include <limits>
#include <cstddef>

#ifndef FLATLOOKUP_HH
#define FLATLOOKUP_HH

#include <numeric/types.hh>
#include <vector>

namespace numeric
{
namespace coordinate_fitting
{

template <class QueryType, class EntryType, class Real=double>
class FlatLookup
{
public:
	~FlatLookup() {};

	template <class InputIterator>
	void initialize(InputIterator first, InputIterator last)
	{
		entries.assign(first, last);
	}

	// @brief Find first matching object.
	bool first_match(QueryType & query, EntryType & entry, Real & distance)
	{
		prepare_for_query(query);

		distance = std::numeric_limits<Real>::max();
		for ( numeric::Size i = 0; i < entries.size(); i++ ) {
			Real test_distance = entry_distance(query, entries[i]);

			if ( test_distance < distance ) {
				distance = test_distance;
			}
			if ( test_distance < entry_radius(entries[i]) ) {
				distance = test_distance; //redundant.
				entry = entries[i];
				return true;
			}
		}
		return false;
	}

	// @brief Find closest matching object.
	bool closest_match(QueryType & query, EntryType & entry, Real & distance)
	{
		prepare_for_query(query);

		distance = std::numeric_limits<Real>::max();

		for ( numeric::Size i = 0; i < entries.size(); i++ ) {
			Real test_distance = entry_distance(query, entries[i]);
			if (test_distance < distance) { //previous only returned the distance if it was < the entry_radius. This didn't make sense because you want the closest match
				distance = test_distance;
				entry = entries[i];
			}
		}
		if ( distance < std::numeric_limits<Real>::max() ) {
			return true;
		} else {
			return false;
		}
	}


	// @brief Find closest matching object.
	bool closest_match_subset(QueryType & query, EntryType & entry, Real & distance, std::vector<bool> subset)
	{
		prepare_for_query(query);

		distance = std::numeric_limits<Real>::max();

		for ( numeric::Size i = 0; i < entries.size(); i++ ) {
			if(subset[i]){
				Real test_distance = entry_distance(query, entries[i]);
				if (test_distance < distance) { //previous only returned the distance if it was < the entry_radius. This didn't make sense because you want the closest match
					distance = test_distance;
					entry = entries[i];
				}
			}
		}
		if ( distance < std::numeric_limits<Real>::max() ) {
			return true;
		} else {
			return false;
		}
	}


	void all_matches_below_threshold(QueryType & query, std::vector<Size> & entry_loc_v, std::vector<Real> & distance , Real rms_threshold)
	{
		prepare_for_query(query);
		for ( numeric::Size ii = 0; ii < entries.size(); ii++ ) {
			Real test_distance = entry_distance(query, entries[ii]);
			if ( test_distance < rms_threshold ) {
				entry_loc_v.push_back(ii);
				distance.push_back(test_distance);
			}
		}
	}


	// @brief Called before query evaluation.
	virtual void prepare_for_query(QueryType &) { }

	// @brief Calculate distance between query object and entry
	virtual Real entry_distance(QueryType & q, EntryType & e2) = 0;

	// @brief Get maximum radius from entry
	virtual Real entry_radius(EntryType & e) = 0;

	// @brief Souce entries query
	std::vector<EntryType>  entries;
};

}
}

#endif
