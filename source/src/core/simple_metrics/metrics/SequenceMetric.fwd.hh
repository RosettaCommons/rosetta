// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/SequenceMetric.fwd.hh
/// @brief A SimpleMetric to output the single-letter OR three-letter sequence of a protein or subset of positions/regions using a ResidueSelector.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org) -- Added support for writing full residue type names or basenames.

#ifndef INCLUDED_core_simple_metrics_metrics_SequenceMetric_fwd_hh
#define INCLUDED_core_simple_metrics_metrics_SequenceMetric_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace simple_metrics {
namespace metrics {

/// @brief The mode for this metric.  If you add to this list, be sure to update the map associating this enum
/// with corresponding strings in the .cc file.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
enum class SequenceMetricMode {
	ONELETTER_CODE=1,
	THREELETTER_CODE,
	BASE_NAME,
	FULL_NAME,
	INVALID_MODE, //Keep this second-to-last.
	END_OF_LIST = INVALID_MODE //Keep this last
};

class SequenceMetric;

typedef utility::pointer::shared_ptr< SequenceMetric > SequenceMetricOP;
typedef utility::pointer::shared_ptr< SequenceMetric const > SequenceMetricCOP;

} //core
} //simple_metrics
} //metrics

#endif //INCLUDED_core_simple_metrics_metrics_SequenceMetric_fwd_hh
