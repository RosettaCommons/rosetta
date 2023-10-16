// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/composite_metrics/BestMutationsFromProbabilitiesMetric.fwd.hh
/// @brief A class for calculating the mutations with the highest delta_probability to the current residues from a PerResidueProbabilitiesMetric.
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)

#ifndef INCLUDED_core_simple_metrics_composite_metrics_BestMutationsFromProbabilitiesMetric_fwd_hh
#define INCLUDED_core_simple_metrics_composite_metrics_BestMutationsFromProbabilitiesMetric_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Forward
namespace core {
namespace simple_metrics {
namespace composite_metrics {

class BestMutationsFromProbabilitiesMetric;

using BestMutationsFromProbabilitiesMetricOP = utility::pointer::shared_ptr< BestMutationsFromProbabilitiesMetric >;
using BestMutationsFromProbabilitiesMetricCOP = utility::pointer::shared_ptr< BestMutationsFromProbabilitiesMetric const >;

} //composite_metrics
} //simple_metrics
} //core

#endif //INCLUDED_core_simple_metrics_composite_metrics_BestMutationsFromProbabilitiesMetric_fwd_hh
