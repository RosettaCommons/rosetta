// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/FilterValueMetric.fwd.hh
/// @brief Convert the result of a Filter's report_sm() to a SimpleMetric
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_filters_FilterValueMetric_fwd_hh
#define INCLUDED_protocols_filters_FilterValueMetric_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace filters {

class FilterValueMetric;

using FilterValueMetricOP = utility::pointer::shared_ptr< FilterValueMetric >;
using FilterValueMetricCOP = utility::pointer::shared_ptr< FilterValueMetric const >;

} //filters
} //protocols

#endif //INCLUDED_protocols_filters_FilterValueMetric_fwd_hh
