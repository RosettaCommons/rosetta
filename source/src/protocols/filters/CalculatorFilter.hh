// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/CalculatorFilter.hh
/// @brief Combine several filters in a (semi) arbitrary calculation
/// @author Rocco Moretti (rmoretti@u.washington.edu)

#ifndef INCLUDED_protocols_filters_CalculatorFilter_hh
#define INCLUDED_protocols_filters_CalculatorFilter_hh

// Unit Headers
#include <protocols/filters/CalculatorFilter.fwd.hh>

// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/simple_metrics/metrics/CalculatorMetric.fwd.hh>

#include <string>

// Unit headers

namespace protocols {
namespace filters {

class CalculatorFilter : public protocols::filters::Filter
{
public:
	CalculatorFilter();
	CalculatorFilter(std::string const & equation);

	bool apply( core::pose::Pose const & pose ) const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Real compute( core::pose::Pose const & pose ) const;

	protocols::filters::FilterOP clone() const override {
		return utility::pointer::make_shared< CalculatorFilter >( *this );
	}
	protocols::filters::FilterOP fresh_instance() const override{
		return utility::pointer::make_shared< CalculatorFilter >();
	}

	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &
	) override;

	void threshold( core::Real threshold) { threshold_ = threshold; }

	void add_filter( std::string const & name, protocols::filters::FilterOP filter );

	void add_reported_value( std::string const & name, std::string const & report_key );

	void add_constant( std::string const & name, core::Real value );

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	core::simple_metrics::metrics::CalculatorMetricOP metric_;
	core::Real threshold_;

};

} // filters
} // protocols

#endif //INCLUDED_protocols_Filters_CalculatorFilter_HH_
