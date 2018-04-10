// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/SasaMetricCreator.hh
/// @brief Simple metrics for calculating and adding to pose
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_metrics_test_classes_creators_hh
#define INCLUDED_core_simple_metrics_metrics_test_classes_creators_hh

// Unit headers
#include <core/simple_metrics/SimpleMetricCreator.hh>

// Protocol headers
#include <core/simple_metrics/SimpleMetric.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace core {
namespace simple_metrics {

class TestStringMetricCreator : public core::simple_metrics::SimpleMetricCreator {
public:


	/// @brief Instantiate a particular SimpleMetric
	virtual SimpleMetricOP
	create_simple_metric() const override;

	/// @brief Return a string that will be used to instantiate the particular SimpleMetric
	virtual std::string
	keyname() const override;

	virtual void
	provide_xml_schema( utility::tag::XMLSchemaDefinition &) const override;
};

class TestIntegerMetricCreator : public core::simple_metrics::SimpleMetricCreator {
public:


	/// @brief Instantiate a particular SimpleMetric
	virtual SimpleMetricOP
	create_simple_metric() const override;

	/// @brief Return a string that will be used to instantiate the particular SimpleMetric
	virtual std::string
	keyname() const override;

	virtual void
	provide_xml_schema( utility::tag::XMLSchemaDefinition &) const override;
};

class TestRealMetricCreator : public core::simple_metrics::SimpleMetricCreator {
public:


	/// @brief Instantiate a particular SimpleMetric
	virtual SimpleMetricOP
	create_simple_metric() const override;

	/// @brief Return a string that will be used to instantiate the particular SimpleMetric
	virtual std::string
	keyname() const override;

	virtual void
	provide_xml_schema( utility::tag::XMLSchemaDefinition &) const override;
};

class TestCompositeStringMetricCreator : public core::simple_metrics::SimpleMetricCreator {
public:


	/// @brief Instantiate a particular SimpleMetric
	virtual SimpleMetricOP
	create_simple_metric() const override;

	/// @brief Return a string that will be used to instantiate the particular SimpleMetric
	virtual std::string
	keyname() const override;

	virtual void
	provide_xml_schema( utility::tag::XMLSchemaDefinition &) const override;
};

class TestCompositeIntegerMetricCreator : public core::simple_metrics::SimpleMetricCreator {
public:


	/// @brief Instantiate a particular SimpleMetric
	virtual SimpleMetricOP
	create_simple_metric() const override;

	/// @brief Return a string that will be used to instantiate the particular SimpleMetric
	virtual std::string
	keyname() const override;

	virtual void
	provide_xml_schema( utility::tag::XMLSchemaDefinition &) const override;
};

class TestCompositeRealMetricCreator : public core::simple_metrics::SimpleMetricCreator {
public:


	/// @brief Instantiate a particular SimpleMetric
	virtual SimpleMetricOP
	create_simple_metric() const override;

	/// @brief Return a string that will be used to instantiate the particular SimpleMetric
	virtual std::string
	keyname() const override;

	virtual void
	provide_xml_schema( utility::tag::XMLSchemaDefinition &) const override;
};

} //simple_metrics
} //core

#endif //INCLUDED_core_simple_metrics_SasaMetric_fwd_hh
