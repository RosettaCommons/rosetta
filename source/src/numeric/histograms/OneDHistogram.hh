// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//////////////////////////////////////////////////////////////////////
///
/// @brief
/// A 1D histogram based upon a map structure
///
/// @details
/// Very simple class for histograms based upon maps. You provide the key, which is templated,
/// meaning that the key can be a string, real, size, enum. It will return a count, if you want it
///
///
///
/// @author Steven Combs
///
/////////////////////////////////////////////////////////////////////////


#ifndef INCLUDED_numeric_histograms_OneDHistogram_hh
#define INCLUDED_numeric_histograms_OneDHistogram_hh


#include <map>
#include <platform/types.hh>

namespace numeric{
namespace histograms{


template<typename key1>
class OneDHistogram {

public:

	OneDHistogram<key1>()= default;

	void insert_data(key1 key_1, platform::Size counts){
		histogram_.insert(std::make_pair(key_1, counts));
	}

	platform::Size lookup_counts(key1 key_1){
		platform::Size counts(histogram_.find(key_1)->second );

		return counts;
	}


private:
std::map< key1, platform::Size > histogram_;


};


}
}
#endif /* INCLUDED_numeric_histograms_OneDHistogram_hh */
