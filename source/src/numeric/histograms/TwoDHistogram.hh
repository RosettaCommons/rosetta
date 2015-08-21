// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
///
/// @brief
/// A 2D histogram based upon a map structure
///
/// @details
/// Very simple class for histograms based upon maps. You provide the key, which is templated,
/// meaning that the two keys can be strings, reals, sizes. It will return a count, if you want it
///
///
///
/// @author Steven Combs
///
/////////////////////////////////////////////////////////////////////////


#ifndef INCLUDED_numeric_histograms_TwoDHistogram_hh
#define INCLUDED_numeric_histograms_TwoDHistogram_hh

#include <map>
#include <platform/types.hh>


namespace numeric {
namespace histograms {


template<typename key1, typename key2>

class TwoDHistogram {


public:
	TwoDHistogram() {

	}

	/// @brief increase the counts of a histogram using paired keys
	void increase_count(std::pair<key1, key2> paired_keys){
		if ( histogram_.find( paired_keys) != histogram_.end() ) {
			++histogram_[ paired_keys];
		} else {
			histogram_.insert(std::make_pair(paired_keys, 1));
		}
	}

	/// @brief overflow function to increase count taking 2 keys to increment counts of histogram
	void increase_count(key1 key_1, key2 key_2){
		std::pair<key1, key2> paired_keys(key_1, key_2);
		if ( histogram_.find( paired_keys) != histogram_.end() ) {
			++histogram_[ paired_keys];
		} else {
			histogram_.insert(std::make_pair<paired_keys, 1>);
		}
	}


	/// @brief insert data to the histogram
	void insert_data(key1 key_1, key2 key_2, platform::Size counts){
		std::pair<key1, key2> paired_keys(key_1, key_2);
		histogram_.insert(std::make_pair(paired_keys, counts));
	}

	/// @brief overload function to insert data to the histogram
	void insert_data(std::pair<key1, key2 > paired_keys, platform::Size counts){
		histogram_.insert(std::make_pair(paired_keys, counts));
	}

	/// @brief look up data based upon your keys
	platform::Size lookup_counts(key1 key_1, key2 key_2){
		std::pair<key1, key2> paired_keys(key_1, key_2);

		platform::Size counts(histogram_.find(paired_keys)->second);

		return counts;
	}

	/// @brief overload function to lookup counts
	platform::Size lookup_counts(std::pair<key1, key2> paired_keys){
		platform::Size counts(histogram_.find(paired_keys)->second);

		return counts;
	}


private:
	std::map< std::pair< key1, key2>, platform::Size > histogram_;


};


}
}


#endif /* INCLUDED_numeric_histograms_TwoDHistogram_hh */
