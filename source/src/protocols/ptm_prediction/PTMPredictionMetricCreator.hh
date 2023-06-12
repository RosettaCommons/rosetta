// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/ptm_prediction/PTMPredictionMetricCreator.hh
/// @brief PTM prediction
/// @author MoritzErtelt (moritz.ertelt@googlemail.com)

#ifndef INCLUDED_protocols_ptm_prediction_PTMPredictionMetricCreator_HH
#define INCLUDED_protocols_ptm_prediction_PTMPredictionMetricCreator_HH

// Unit headers
#include <core/simple_metrics/SimpleMetricCreator.hh>

// Protocol headers
#include <core/simple_metrics/SimpleMetric.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace ptm_prediction {

class PTMPredictionMetricCreator : public core::simple_metrics::SimpleMetricCreator {
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

} //ptm_prediction
} //protocols

#endif //INCLUDED_protocols_ptm_prediction_PTMPredictionMetricCreator_HH

