// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/analysis/simple_metrics/ConstraintsMetricCreator.hh
/// @brief A simple metric that writes out all the constraints in a pose or sub-region of a pose,
/// in a format matching the (non-enzdes) constraints file format.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_analysis_simple_metrics_ConstraintsMetricCreator_HH
#define INCLUDED_protocols_analysis_simple_metrics_ConstraintsMetricCreator_HH

// Unit headers
#include <core/simple_metrics/SimpleMetricCreator.hh>

// Protocol headers
#include <core/simple_metrics/SimpleMetric.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace analysis {
namespace simple_metrics {

class ConstraintsMetricCreator : public core::simple_metrics::SimpleMetricCreator {
public:


	/// @brief Instantiate a particular SimpleMetric
	core::simple_metrics::SimpleMetricOP
	create_simple_metric() const override;

	/// @brief Return a string that will be used to instantiate the particular SimpleMetric
	std::string
	keyname() const override;

	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition &) const override;
};

} //simple_metrics
} //analysis
} //protocols

#endif //INCLUDED_protocols_analysis_simple_metrics_ConstraintsMetricCreator_HH

