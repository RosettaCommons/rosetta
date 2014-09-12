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

/// @brief Declaration of a filter to exclude a given set of residues
/// @author Andrea Bazzoli (bazzoli@ku.edu)

#ifndef INCLUDED_ExcludedFilter_hh
#define INCLUDED_ExcludedFilter_hh

#include "ResidueMask.hh"

namespace devel {
namespace constel {

class ExcludedFilter {

	// excluded[i] is true iff the ith residue in the pose must be excluded from any
	// constellation (i=1,...,N, where N is the number of residues in the pose)
	static ResidueMaskOP excluded;

	public:
	static bool hasnt_excluded(Pose const&, utility::vector1<Size> const& cnl);
	static void init(Pose& ps, std::string const& ex_fname);
};

} // constel
} // devel

#endif
