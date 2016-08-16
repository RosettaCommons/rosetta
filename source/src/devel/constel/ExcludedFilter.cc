// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief Definition of class ExcludedFilter
/// @author Andrea Bazzoli (bazzoli@ku.edu)

#include "ExcludedFilter.hh"

namespace devel {
namespace constel {


// allocation of the mask of excluded residues
ResidueMaskOP ExcludedFilter::excluded;


/// @brief: reads from file the mask of excluded residues
///
/// @param[in] ps pose for whose residues the mask is built
/// @param[in] ex_fname path to the file describing the mask
///
void ExcludedFilter::init(Pose& ps, std::string const& ex_fname) {

	excluded = ResidueMaskOP( new ResidueMask(ps, ex_fname) );
}


/// @brief: given a constellation, returns true if it doesn't contain any
///  residue to be excluded from any constellation; returns false otherwise.
///
/// @param[in] ps pose to which the constellation belongs (dummy)
/// @param[in] cnl indexes in the pose of the residues forming the
///  constellation
///
bool ExcludedFilter::hasnt_excluded(Pose const&,
	utility::vector1<Size> const& cnl) {

	Size const N = cnl.size();
	for ( Size i=1; i<=N; ++i ) {
		if ( (*excluded)[cnl[i]] ) {
			return false;
		}
	}

	return true;
}

} // constel
} // devel
