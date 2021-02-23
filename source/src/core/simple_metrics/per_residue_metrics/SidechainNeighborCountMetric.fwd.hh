// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/SidechainNeighborCountMetric.fwd.hh
/// @brief A metric for calculating each sidechains neighbors based on cones.  This metric uses the same core code as the LayerSelector.  It can be combined with the SimpleMetricSelector to select any set of residues based on burial.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_per_residue_metrics_SidechainNeighborCountMetric_fwd_hh
#define INCLUDED_core_simple_metrics_per_residue_metrics_SidechainNeighborCountMetric_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace simple_metrics {
namespace per_residue_metrics {

class SidechainNeighborCountMetric;

using SidechainNeighborCountMetricOP = utility::pointer::shared_ptr< SidechainNeighborCountMetric >;
using SidechainNeighborCountMetricCOP = utility::pointer::shared_ptr< SidechainNeighborCountMetric const >;

} //per_residue_metrics
} //simple_metrics
} //core

#endif //INCLUDED_core_simple_metrics_per_residue_metrics_SidechainNeighborCountMetric_fwd_hh
