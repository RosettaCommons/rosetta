// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/DihedralDistanceMetric.fwd.hh
/// @brief A metric to calculate the dihedral distance between two poses or the input and the set cmd-line native.  Can set a subset of residues to calculate via ResidueSelector.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_metrics_DihedralDistanceMetric_fwd_hh
#define INCLUDED_core_simple_metrics_metrics_DihedralDistanceMetric_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace simple_metrics {
namespace metrics {

class DihedralDistanceMetric;

typedef utility::pointer::shared_ptr< DihedralDistanceMetric > DihedralDistanceMetricOP;
typedef utility::pointer::shared_ptr< DihedralDistanceMetric const > DihedralDistanceMetricCOP;

} //core
} //simple_metrics
} //metrics

#endif //INCLUDED_protocols_analysis_simple_metrics_DihedralDistanceMetric_fwd_hh
