// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/saxs/DistanceHistogram.fwd.hh
/// @brief Forward declaration of a DistanceHistogram data type
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_core_scoring_saxs_DistanceHistogram_fwd_hh
#define INCLUDED_core_scoring_saxs_DistanceHistogram_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace saxs {

class DistanceHistogram;

typedef utility::pointer::shared_ptr<DistanceHistogram>
	DistanceHistogramOP;
typedef utility::pointer::shared_ptr<DistanceHistogram const>
	DistanceHistogramCOP;

} // saxs
} // scoring
} // core


#endif /* INCLUDED_core_scoring_saxs_DistanceHistogram_FWD_HH */
