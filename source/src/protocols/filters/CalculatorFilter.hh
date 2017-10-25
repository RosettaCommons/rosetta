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
#include <protocols/moves/Mover.fwd.hh>

#include <numeric/Calculator.fwd.hh>
#include <utility/vector1.hh>

#include <string>

// Unit headers

namespace protocols {
namespace filters {

class CalculatorFilter : public protocols::filters::Filter
{
public:
	CalculatorFilter();
	CalculatorFilter(std::string equation);
	CalculatorFilter(CalculatorFilter const & other);
	~CalculatorFilter() override;

	bool apply( core::pose::Pose const & pose ) const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Real compute( core::pose::Pose const & pose ) const;

	protocols::filters::FilterOP clone() const override {
		return protocols::filters::FilterOP( new CalculatorFilter( *this ) );
	}
	protocols::filters::FilterOP fresh_instance() const override{
		return protocols::filters::FilterOP( new CalculatorFilter() );
	}

	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & ) override;
	void threshold( core::Real threshold) { threshold_ = threshold; }

	void add_filter( std::string name, protocols::filters::FilterOP filter );

	void add_reported_value( std::string name, std::string report_key );

	void add_constant( std::string name, core::Real value );

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	numeric::CalculatorOP calc_;
	std::map<std::string, core::Real> values_;
	std::map<std::string, protocols::filters::FilterOP> filters_;
	std::map<std::string, std::string> reported_values_;
	core::Real threshold_;

};

} // filters
} // protocols

#endif //INCLUDED_protocols_Filters_CalculatorFilter_HH_
