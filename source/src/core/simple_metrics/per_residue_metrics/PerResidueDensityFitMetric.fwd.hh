// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/PerResidueDensityFitMetric.fwd.hh
/// @brief A per-residue metric that will calculate the density fit for each residue using either a correlation or a zscore.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_per_residue_metrics_PerResidueDensityFitMetric_fwd_hh
#define INCLUDED_core_simple_metrics_per_residue_metrics_PerResidueDensityFitMetric_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace simple_metrics {
namespace per_residue_metrics {

class PerResidueDensityFitMetric;

typedef utility::pointer::shared_ptr< PerResidueDensityFitMetric > PerResidueDensityFitMetricOP;
typedef utility::pointer::shared_ptr< PerResidueDensityFitMetric const > PerResidueDensityFitMetricCOP;

} //core
} //simple_metrics
} //per_residue_metrics

#endif //INCLUDED_core_simple_metrics_per_residue_metrics_PerResidueDensityFitMetric_fwd_hh
