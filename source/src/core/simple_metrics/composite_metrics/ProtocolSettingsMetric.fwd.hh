// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/composite_metrics/ProtocolSettingsMetric.fwd.hh
/// @brief This Metric reports options that have been set in the command line and splits script_vars.  Each option name is the type and the setting is the value in the map.  This is primarily aimed at benchmarking and record-keeping for large-scale rosetta runs or experiments.  It works with both the global and local OptionsCollection to enable its use in JD3.  It is analogous to the ProtocolFeatures reporter, with more options for xml-based variables.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_composite_metrics_ProtocolSettingsMetric_fwd_hh
#define INCLUDED_core_simple_metrics_composite_metrics_ProtocolSettingsMetric_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace simple_metrics {
namespace composite_metrics {

class ProtocolSettingsMetric;

typedef utility::pointer::shared_ptr< ProtocolSettingsMetric > ProtocolSettingsMetricOP;
typedef utility::pointer::shared_ptr< ProtocolSettingsMetric const > ProtocolSettingsMetricCOP;

} //core
} //simple_metrics
} //composite_metrics

#endif //INCLUDED_core_simple_metrics_composite_metrics_ProtocolSettingsMetric_fwd_hh
