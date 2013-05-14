// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @brief Definition of class MasterFilter
/// @author Andrea Bazzoli

#include <devel/constel/MasterFilter.hh>


namespace devel {
namespace constel {

/// @brief Array of filters applied to a constellation.
utility::vector1<MasterFilter::FiltPtr> MasterFilter::filters;

///
/// @brief Adds a filter to the array of filters.
///
/// @param[in] f pointer to the filter.
///
void MasterFilter::addfilt(FiltPtr f) {

  filters.push_back(f);
}


///
/// @brief Tells whether a constellation is valid.
///
/// @param[in] ps pose to which all residues forming the constellation belong.
/// @param[in] cnl indexes in the pose of the residues forming the
/// 	constellation.
///
/// @return true if the constellation is valid, false otherwise.
///
/// @details Currently, a constellation is deemed to be valid if and only if it
/// 	passes all the filters. It is up to the user to decide which filters are
/// 	to be applied.
///
bool MasterFilter::is_constel_valid(Pose const& ps,
	utility::vector1<Size> const& cnl) {

  for(unsigned int i=1; i<=filters.size(); ++i)
    if( !filters[i](ps, cnl) )
      return false;

  return true;
}

} // constel
} // devel
