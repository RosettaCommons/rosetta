// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
