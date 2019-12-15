// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constel/ExcludedFilter.hh
/// @brief Declaration of a filter to exclude a given set of residues
/// @author Andrea Bazzoli

#ifndef INCLUDED_ExcludedFilter_hh
#define INCLUDED_ExcludedFilter_hh

#include "ResidueMask.hh"

namespace protocols {
namespace constel {

class ExcludedFilter {

	// excluded[i] is true iff the ith residue in the pose must be excluded from any
	// constellation (i=1,...,N, where N is the number of residues in the pose)
	static ResidueMaskOP excluded;

public:
	static bool hasnt_excluded(core::pose::Pose const&, utility::vector1<core::Size> const& cnl);
	static void init(core::pose::Pose& ps, std::string const& ex_fname);
};

} // constel
} // protocols

#endif
