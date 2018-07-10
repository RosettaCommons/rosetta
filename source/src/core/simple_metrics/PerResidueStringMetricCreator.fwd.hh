// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/simple_metrics/PerResidueStringMetricCreator.fwd.hh
/// @brief  Forward declaration of a class that instantiates a particular PerResidueStringMetric
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_PerResidueStringMetricCreator_FWD_HH
#define INCLUDED_core_simple_metrics_PerResidueStringMetricCreator_FWD_HH

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace simple_metrics {

class PerResidueStringMetricCreator;

typedef utility::pointer::shared_ptr< PerResidueStringMetricCreator > PerResidueStringMetricCreatorOP;
typedef utility::pointer::shared_ptr< PerResidueStringMetricCreator const > PerResidueStringMetricCreatorCOP;

} //namespace simple_metrics
} //namespace core

#endif
