// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/esm_perplexity/PseudoPerplexityMetricCreator.hh
/// @brief Predicting amino acid probabilities using the ESM language model
/// @author MoritzErtelt (moritz.ertelt@googlemail.com)

#ifndef INCLUDED_protocols_esm_perplexity_EsmPerplexityMetricCreator_HH
#define INCLUDED_protocols_esm_perplexity_EsmPerplexityMetricCreator_HH

// Unit headers
#include <core/simple_metrics/SimpleMetricCreator.hh>

// Protocol headers
#include <core/simple_metrics/SimpleMetric.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace esm_perplexity {

class PseudoPerplexityMetricCreator : public core::simple_metrics::SimpleMetricCreator {
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

} //esm_perplexity
} //protocols

#endif //INCLUDED_protocols_esm_perplexity_EsmPerplexityMetricCreator_HH

