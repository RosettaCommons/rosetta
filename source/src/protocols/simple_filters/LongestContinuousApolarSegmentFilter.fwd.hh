// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/LongestContinuousApolarSegmentFilter.fwd.hh
/// @brief This filter computes the longest continuous stretch of polar residues within a pose or selection.
/// @author Yang Hsia (yhsia@uw.edu)

#ifndef INCLUDED_protocols_simple_filters_LongestContinuousApolarSegmentFilter_fwd_hh
#define INCLUDED_protocols_simple_filters_LongestContinuousApolarSegmentFilter_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace simple_filters {

class LongestContinuousApolarSegmentFilter;

typedef utility::pointer::shared_ptr< LongestContinuousApolarSegmentFilter > LongestContinuousApolarSegmentFilterOP;
typedef utility::pointer::shared_ptr< LongestContinuousApolarSegmentFilter const > LongestContinuousApolarSegmentFilterCOP;

} //protocols
} //simple_filters

#endif //INCLUDED_protocols_simple_filters_LongestContinuousApolarSegmentFilter_fwd_hh
