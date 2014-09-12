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

/// @brief Definition of class ExcludedFilter
/// @author Andrea Bazzoli (bazzoli@ku.edu)

#include "ExcludedFilter.hh"

namespace devel {
namespace constel {


// allocation of the mask of excluded residues
ResidueMaskOP ExcludedFilter::excluded;


///
/// @brief: reads from file the mask of excluded residues
///
/// @param[in] ps pose for whose residues the mask is built
/// @param[in] ex_fname path to the file describing the mask
///
void ExcludedFilter::init(Pose& ps, std::string const& ex_fname) {

	excluded = new ResidueMask(ps, ex_fname);
}


///
/// @brief: given a constellation, returns true if it doesn't contain any
/// 	residue to be excluded from any constellation; returns false otherwise.
///
/// @param[in] ps pose to which the constellation belongs (dummy)
/// @param[in] cnl indexes in the pose of the residues forming the
/// 	constellation
///
bool ExcludedFilter::hasnt_excluded(Pose const&,
	utility::vector1<Size> const& cnl) {

	Size const N = cnl.size();
	for(Size i=1; i<=N; ++i)
		if((*excluded)[cnl[i]])
			return false;

	return true;
}

} // constel
} // devel
