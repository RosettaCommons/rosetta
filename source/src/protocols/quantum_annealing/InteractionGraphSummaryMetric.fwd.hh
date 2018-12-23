// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/quantum_annealing/InteractionGraphSummaryMetric.fwd.hh
/// @brief A simple metric that allows an InteractionGraph to be written out in a format that external annealers (including quantum annealers) can read.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_quantum_annealing_InteractionGraphSummaryMetric_fwd_hh
#define INCLUDED_protocols_quantum_annealing_InteractionGraphSummaryMetric_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace quantum_annealing {

class InteractionGraphSummaryMetric;

typedef utility::pointer::shared_ptr< InteractionGraphSummaryMetric > InteractionGraphSummaryMetricOP;
typedef utility::pointer::shared_ptr< InteractionGraphSummaryMetric const > InteractionGraphSummaryMetricCOP;

} //protocols
} //quantum_annealing

#endif //INCLUDED_protocols_quantum_annealing_InteractionGraphSummaryMetric_fwd_hh
