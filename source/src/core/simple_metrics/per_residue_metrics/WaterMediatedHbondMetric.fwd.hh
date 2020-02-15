// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/WaterMediatedHbondMetric.fwd.hh
/// @brief A metric to measure hydrogen bonds between a set of residues that are water-mediated.  Depth of 1 is default where one water mediates the interface between these residues.  Depth can be set.  Make sure to use the -include_waters flag to have Rosetta not ignore HOH

/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_per_residue_metrics_WaterMediatedHbondMetric_fwd_hh
#define INCLUDED_core_simple_metrics_per_residue_metrics_WaterMediatedHbondMetric_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace simple_metrics {
namespace per_residue_metrics {

class WaterMediatedHbondMetric;

using WaterMediatedHbondMetricOP = utility::pointer::shared_ptr< WaterMediatedHbondMetric >;
using WaterMediatedHbondMetricCOP = utility::pointer::shared_ptr< WaterMediatedHbondMetric const >;

} //per_residue_metrics
} //simple_metrics
} //core

#endif //INCLUDED_core_simple_metrics_per_residue_metrics_WaterMediatedHbondMetric_fwd_hh
