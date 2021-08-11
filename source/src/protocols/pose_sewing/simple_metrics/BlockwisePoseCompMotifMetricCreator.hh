// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/simple_metrics/BlockwisePoseCompMotifMetricCreator.hh
/// @brief composite metric returning motif scores that pass per block
/// @author frankdt (frankdt@email.unc.edu)

#ifndef INCLUDED_protocols_pose_sewing_simple_metrics_BlockwisePoseCompMotifMetricCreator_HH
#define INCLUDED_protocols_pose_sewing_simple_metrics_BlockwisePoseCompMotifMetricCreator_HH

// Unit headers
#include <core/simple_metrics/SimpleMetricCreator.hh>

// Protocol headers
#include <core/simple_metrics/SimpleMetric.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace pose_sewing {
namespace simple_metrics {

class BlockwisePoseCompMotifMetricCreator : public core::simple_metrics::SimpleMetricCreator {
public:


	/// @brief Instantiate a particular SimpleMetric
	virtual core::simple_metrics::SimpleMetricOP
	create_simple_metric() const override;

	/// @brief Return a string that will be used to instantiate the particular SimpleMetric
	virtual std::string
	keyname() const override;

	virtual void
	provide_xml_schema( utility::tag::XMLSchemaDefinition &) const override;
};

} //simple_metrics
} //pose_sewing
} //protocols

#endif //INCLUDED_protocols_pose_sewing_simple_metrics_BlockwisePoseCompMotifMetricCreator_HH

