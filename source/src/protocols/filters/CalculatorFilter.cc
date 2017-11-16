// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/CalculatorFilter.cc
/// @brief Combine several filters in a (semi) arbitrary calculation
/// @author Rocco Moretti (rmoretti@u.washington.edu)

#include <protocols/filters/CalculatorFilter.hh>
#include <protocols/filters/CalculatorFilterCreator.hh>

#include <core/pose/Pose.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <utility/vector1.hh>
#include <numeric/Calculator.hh>
#include <numeric/random/random.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


namespace protocols {
namespace filters {

static basic::Tracer TR( "protocols.filters.CalculatorFilter" );

/// @brief default ctor
CalculatorFilter::CalculatorFilter() :
	Filter( "CalculatorFilter" ),
	threshold_(0)
{}

CalculatorFilter::CalculatorFilter(std::string equation) :
	Filter( "CalculatorFilter" ),
	threshold_(0)
{
	calc_ = numeric::CalculatorOP( new numeric::Calculator(equation) );
}

CalculatorFilter::CalculatorFilter(CalculatorFilter const & other) :
	Filter( "CalculatorFilter" ),
	calc_(other.calc_),
	values_(other.values_),
	filters_(other.filters_),
	threshold_(other.threshold_)
{}


CalculatorFilter::~CalculatorFilter() = default;

bool
CalculatorFilter::apply(core::pose::Pose const & pose) const
{
	core::Real value = compute(pose);
	return value <= threshold_;
}

core::Real
CalculatorFilter::report_sm( core::pose::Pose const & pose) const
{
	return( compute(pose) );
}

void
CalculatorFilter::report( std::ostream & out, core::pose::Pose const & pose) const
{
	out<<"CalculatorFilter returns "<< compute(pose) <<std::endl;
}

void
CalculatorFilter::add_filter( std::string name, protocols::filters::FilterOP filter ) {
	if ( ! filter ) {
		utility_exit_with_message("Calculator filter can't use non-existant (null pointer) filter.");
	}
	filters_[name] = filter;
}

void
CalculatorFilter::add_reported_value( std::string name, std::string report_key ) {
	reported_values_[name] = report_key;
}

void
CalculatorFilter::add_constant( std::string name, core::Real value ) {
	values_[name] = value;
}

core::Real
CalculatorFilter::compute(core::pose::Pose const & pose) const {
	debug_assert(calc_);

	std::map< std::string, core::Real > vars(values_);
	for ( auto & filter : filters_ ) {
		debug_assert(filter.second);
		vars[ filter.first ] = (filter.second)->report_sm( pose );
	}
	for ( auto & reported : reported_values_ ) {
		if ( !getPoseExtraScore(pose, reported.second, vars[reported.first]) ) {
			utility_exit_with_message("CalculatorFilter required reported value not yet present in pose.");
		}
	}

	numeric::Real value(999999);
	if ( calc_->compute(vars, value) ) {
		TR.Error << "Problem calculating equation in CalculatorFilter - resultant value likely garbage." << std::endl;
	}
	return value;
}

void
CalculatorFilter::parse_my_tag( utility::tag::TagCOP tag_ptr,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &)
{
	std::string equation = tag_ptr->getOption< std::string >( "equation" );
	threshold_ = tag_ptr->getOption<core::Real>( "threshold", 0.0 );

	for ( utility::tag::TagCOP sub_tag_ptr : tag_ptr->getTags() ) {
		std::string varname( sub_tag_ptr->getOption<std::string>( "name" ) );

		core::Size num_tags = 0;
		for ( auto & name : {"filter", "filter_name", "reported", "value"} ) {
			if ( sub_tag_ptr->hasOption(name) ) {
				num_tags += 1;
			}
		}

		if ( num_tags != 1 ) {
			utility_exit_with_message("CalculatorFilter subtag must have one of 'filter', 'reported', or 'value' set.");
		}

		if ( sub_tag_ptr->hasOption("filter") || sub_tag_ptr->hasOption("filter_name") ) {
			std::string filter_name;
			if ( sub_tag_ptr->hasOption("filter_name") ) {
				filter_name = sub_tag_ptr->getOption<std::string>( "filter_name" );
			} else {
				filter_name = sub_tag_ptr->getOption<std::string>( "filter" );
			}
			add_filter( varname, protocols::rosetta_scripts::parse_filter( filter_name , filters ) );

			continue;
		} else if ( sub_tag_ptr->hasOption("reported") ) {
			add_reported_value( varname, sub_tag_ptr->getOption<std::string>( "reported" ) );
		} else if ( sub_tag_ptr->hasOption("value") ) {
			add_constant( varname, sub_tag_ptr->getOption<core::Real>( "value" ) );
		} else {
			utility_exit_with_message("CalculatorFilter subtag must have filter, reported, or value set.");
		}
	}

	// Now do a quick sanity check for the equation parsing.
	calc_ = numeric::CalculatorOP( new numeric::Calculator( equation ) );
	std::map< std::string, core::Real > vars(values_);
	for ( auto & filter : filters_ ) {
		vars[ filter.first ] = 1.0 + 0.00001 * numeric::random::uniform(); // Additional random to avoid "1/(alpha - beta)" type situations.
	}
	for ( auto & report : reported_values_ ) {
		vars[ report.first ] = 1.0 + 0.00001 * numeric::random::uniform(); // Additional random to avoid "1/(alpha - beta)" type situations.
	}
	numeric::Real dummy;
	if ( calc_->compute(vars, dummy) ) {
		utility_exit_with_message("Bad equation in CalculatorFilter: " + equation);
	}
}


std::string CalculatorFilter::name() const {
	return class_name();
}

std::string CalculatorFilter::class_name() {
	return "CalculatorFilter";
}

void CalculatorFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	AttributeList subelement_attlist;
	subelement_attlist
		+ XMLSchemaAttribute::required_attribute( "name", xs_string, "Unique name to identify value for use in equation." )
		+ XMLSchemaAttribute( "value", xsct_real, "Specifiy a constant value." )
		+ XMLSchemaAttribute( "filter", xs_string, "Evaluate given filter at calculator evaluation time." )
		+ XMLSchemaAttribute( "filter_name", xs_string, "Evaluate given filter at calculator evaluation time." )
		+ XMLSchemaAttribute( "reported", xs_string, "Retrieve reported value. See 'report_at_end=false' documentation in ParsedProtocol.");

	XMLSchemaSimpleSubelementList subelements;
	subelements
		.add_simple_subelement( "VAR", subelement_attlist, "Specify values to be available in equation." )
		.add_simple_subelement( "Var", subelement_attlist, "Specify values to be available in equation." )
		.add_simple_subelement( "var", subelement_attlist, "Specify values to be available in equation." );

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "equation", xs_string, "Equation to evaluate filter value." )
		+ XMLSchemaAttribute::attribute_w_default( "threshold", xsct_real, "Filter passes if equation value less than threshold, fails otherwise.", "0.0" );

	protocols::filters::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, class_name(), "XRW TO DO", attlist, subelements );
}

std::string CalculatorFilterCreator::keyname() const {
	return CalculatorFilter::class_name();
}

protocols::filters::FilterOP
CalculatorFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new CalculatorFilter );
}

void CalculatorFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CalculatorFilter::provide_xml_schema( xsd );
}


} // filters
} // protocols
