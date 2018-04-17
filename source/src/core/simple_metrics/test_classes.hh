// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/util.hh
/// @brief Test files for any testing of overall SimpleMetric Framework.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_test_classes_hh
#define INCLUDED_core_simple_metrics_test_classes_hh

#include <core/simple_metrics/SimpleMetric.fwd.hh>
#include <core/simple_metrics/RealMetric.hh>
#include <core/simple_metrics/StringMetric.hh>
#include <core/simple_metrics/CompositeRealMetric.hh>
#include <core/simple_metrics/CompositeStringMetric.hh>

#include <core/simple_metrics/test_classes.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Basic headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

//C++ headers
#include <map>

namespace core {
namespace simple_metrics {


///Create subclasses for each type of metric.
// These are used to test downstream functionality, within TestSuites.
class TestStringMetric : public StringMetric {

public:

	TestStringMetric():
		StringMetric(){};

	std::string
	calculate(pose::Pose const &  ) const override{
		return "TESTING";
	};

	std::string
	metric() const override{
		return "SomeString";
	};

	///@brief Name of the class
	std::string
	name() const override{
		return name_static();
	};

	///@brief Name of the class for creator.
	static
	std::string
	name_static();

	SimpleMetricOP
	clone() const override;

public:

	TestStringMetric(TestStringMetric const & ) = default;

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override;

};


class TestRealMetric : public RealMetric {

public:

	TestRealMetric():
		RealMetric(){};

	core::Real
	calculate(pose::Pose const &  ) const override{
		return 1.0;
	};

	std::string
	metric() const override {
		return "SomeReal";
	};

	///@brief Name of the class
	std::string
	name() const override{
		return name_static();
	};

	///@brief Name of the class for creator.
	static
	std::string
	name_static();

	SimpleMetricOP
	clone() const override;

public:

	TestRealMetric(TestRealMetric const & ) = default;

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override;

};


class TestCompositeStringMetric: public CompositeStringMetric {

public:
	TestCompositeStringMetric():
		CompositeStringMetric(){};

	std::map< std::string, std::string>
	calculate(pose::Pose const &  ) const override{
		std::map< std::string, std::string > data;
		data["s_data1"] = "value1";
		data["s_data2"] = "value2";
		return data;
	};

	std::string
	metric() const override {
		return "SomeCompositeString";
	};

	std::string
	name() const override{
		return name_static();
	};

	static
	std::string
	name_static();

	SimpleMetricOP
	clone() const override;

public:

	utility::vector1< std::string >
	get_metric_names() const override{
		utility::vector1< std::string > names;
		names.push_back("s_data1");
		names.push_back("s_data2");
		return names;
	}

	TestCompositeStringMetric(TestCompositeStringMetric const & ) = default;

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override;

};

class TestCompositeRealMetric: public CompositeRealMetric {

public:
	TestCompositeRealMetric():
		CompositeRealMetric(){};

	std::map< std::string, core::Real>
	calculate(pose::Pose const &  ) const override{
		std::map< std::string, core::Real > data;
		data["r_data1"] = 1.0;
		data["r_data2"] = 2.0;
		return data;
	};

	std::string
	metric() const override {
		return "SomeCompositeReal";
	};

	std::string
	name() const override{
		return name_static();
	};

	static
	std::string
	name_static();

	SimpleMetricOP
	clone() const override;

public:

	utility::vector1< std::string >
	get_metric_names() const override{
		utility::vector1< std::string > names;
		names.push_back("r_data1");
		names.push_back("r_data2");
		return names;
	}

	TestCompositeRealMetric(TestCompositeRealMetric const & ) = default;

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override;

};

} //core
} //simple_metrics


#endif //core/simple_metrics_util_hh

