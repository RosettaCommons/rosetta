// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/PerResidueGlycanLayerMetric.fwd.hh
/// @brief A metric that outputs the layer of the glycan tree as measured by the residue distance to the root.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_per_residue_metrics_PerResidueGlycanLayerMetric_fwd_hh
#define INCLUDED_core_simple_metrics_per_residue_metrics_PerResidueGlycanLayerMetric_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace simple_metrics {
namespace per_residue_metrics {

class PerResidueGlycanLayerMetric;

typedef utility::pointer::shared_ptr< PerResidueGlycanLayerMetric > PerResidueGlycanLayerMetricOP;
typedef utility::pointer::shared_ptr< PerResidueGlycanLayerMetric const > PerResidueGlycanLayerMetricCOP;

} //core
} //simple_metrics
} //per_residue_metrics

#endif //INCLUDED_core_simple_metrics_per_residue_metrics_PerResidueGlycanLayerMetric_fwd_hh
