// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/HbondMetric.fwd.hh
/// @brief A metric to report the total h-bonds of a residue, or from a set of residues to another set of residues.  Use the SummaryMetric to get total hbonds of a selection or between selections. See the WaterMediatedBridgedHBondMetric for water-mediated h-bonds.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_per_residue_metrics_HbondMetric_fwd_hh
#define INCLUDED_core_simple_metrics_per_residue_metrics_HbondMetric_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace simple_metrics {
namespace per_residue_metrics {

class HbondMetric;

using HbondMetricOP = utility::pointer::shared_ptr< HbondMetric >;
using HbondMetricCOP = utility::pointer::shared_ptr< HbondMetric const >;

} //per_residue_metrics
} //simple_metrics
} //core

#endif //INCLUDED_core_simple_metrics_per_residue_metrics_HbondMetric_fwd_hh
