// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/SimpleMetricFilter.cc
/// @brief A filter takes any RealMetric and applies a set cutoff to filter the model.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/simple_filters/SimpleMetricFilter.hh>
#include <protocols/simple_filters/SimpleMetricFilterCreator.hh>

#include <core/simple_metrics/SimpleMetricFactory.hh>
#include <core/simple_metrics/SimpleMetric.hh>
#include <core/simple_metrics/RealMetric.hh>
#include <core/simple_metrics/StringMetric.hh>
#include <core/simple_metrics/IntegerMetric.hh>
#include <core/simple_metrics/CompositeRealMetric.hh>
#include <core/simple_metrics/CompositeStringMetric.hh>
#include <core/simple_metrics/CompositeIntegerMetric.hh>
#include <core/simple_metrics/util.hh>
#include <core/pose/Pose.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map

#include <numeric/util.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/util.hh>
#include <utility/string_util.hh>
#include <protocols/filters/filter_schemas.hh>


static basic::Tracer TR( "protocols.simple_filters.SimpleMetricFilter" );

namespace protocols {
namespace simple_filters {

using namespace core::simple_metrics;

SimpleMetricFilter::SimpleMetricFilter():
	protocols::filters::Filter( "SimpleMetricFilter" )
{

}

SimpleMetricFilter::SimpleMetricFilter( core::simple_metrics::SimpleMetricCOP metric, comparison_type co_type):
	protocols::filters::Filter( "SimpleMetricFilter")
{
	metric_ = metric;
	set_comparison_type( co_type );
}

SimpleMetricFilter::SimpleMetricFilter( core::simple_metrics::SimpleMetricCOP metric, comparison_type co_type, std::string composite_action ):
	protocols::filters::Filter( "SimpleMetricFilter")
{
	metric_ = metric;
	set_comparison_type( co_type );
	set_composite_action( composite_action );
}


SimpleMetricFilter::~SimpleMetricFilter()
{}


protocols::filters::FilterOP
SimpleMetricFilter::clone() const
{
	return protocols::filters::FilterOP( new SimpleMetricFilter( *this ) );
}


protocols::filters::FilterOP
SimpleMetricFilter::fresh_instance() const
{
	return protocols::filters::FilterOP( new SimpleMetricFilter );
}

core::Real
SimpleMetricFilter::report_sm( core::pose::Pose const & pose ) const
{
	if ( metric_->type() != "RealMetric" ) {
		TR << "SimpleMetricFilter can only use report_sm for RealMetrics. Returning 0" << std::endl;
		return 0;
	}

	RealMetric const & r_metric = dynamic_cast<RealMetric const & >( *metric_ );

	return r_metric.calculate(pose);
}

///@brief Sets the cutoff type - aka eq, ne, etc. Options are:
///  eq, ne, lt, gt, lt_or_eq, gt_or_eq
void
SimpleMetricFilter::set_comparison_type( comparison_type co_type ){
	type_ = co_type;
}

void
SimpleMetricFilter::set_composite_action(std::string const & composite_action ){
	composite_action_ = composite_action;
}

///@brief Set the SimpleMetric (RealMetric) that we will be using to filter.
void
SimpleMetricFilter::set_simple_metric( core::simple_metrics::SimpleMetricCOP metric ){
	metric_ = metric;
}

///@brief Set the cutoff value.
void
SimpleMetricFilter::set_cutoff( core::Real cutoff ){
	cutoff_ = cutoff;
}

void
SimpleMetricFilter::set_match_string(std::string const & match_string ){
	match_ = match_string;
}

///@brief Set the sigfigs we will use in our comparisons.
/// Default is .0001
void
SimpleMetricFilter::set_epsilon( core::Real epsilon){
	epsilon_ = epsilon;
}


std::string SimpleMetricFilter::name() const {
	return class_name();
}

std::string SimpleMetricFilter::class_name() {
	return "SimpleMetricFilter";
}

std::map< std::string, comparison_type >
get_string_comparison_type_map(){

	std::map< std::string, comparison_type > types;
	types["eq"] = eq;
	types["ne"] = ne;
	types["lt"] = lt;
	types["gt"] = gt;
	types["lt_or_eq"] = lt_or_eq;
	types["gt_or_eq"] = gt_or_eq;

	return types;
}

utility::vector1< std::string >
get_string_comparison_type_strings(){
	std::map< std::string, comparison_type > types = get_string_comparison_type_map();
	utility::vector1< std::string > type_list;
	for ( auto type_pair : types ) {
		type_list.push_back(type_pair.first);
	}
	return type_list;
}

void
SimpleMetricFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	//We should also accept a comma-separated list of previously defined simple metrics
	// Full name of simple_metrics as this is exactly what they are.

	SimpleMetricCOP metric = get_metric_from_datamap_and_subtags(tag, data);
	set_simple_metric(metric);

	std::map< std::string, comparison_type > types = get_string_comparison_type_map();
	set_comparison_type( types[tag->getOption< std::string >("comparison_type")]);
	set_epsilon( tag->getOption< core::Real >("epsilon", epsilon_));
	set_cutoff( tag->getOption< core::Real >("cutoff", 0));

	if ( tag->hasOption("composite_action") ) {
		set_composite_action( tag->getOption< std::string >("composite_action"));
	}
	if ( tag->hasOption("match") ) {
		set_match_string( tag->getOption< std::string >("match"));
	}

}

void SimpleMetricFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute::attribute_w_default("epsilon", xsct_real , "Epsilon for numerical comparisons", ".0001");
	attlist + XMLSchemaAttribute( "metric", xs_string, "The metric to run in this filter.  Currently must be a RealMetric. IE it calculates one float value.  ");

	utility::vector1< std::string > types = get_string_comparison_type_strings();
	utility::tag::add_schema_restrictions_for_strings( xsd, "comparison_types", types);


	attlist + XMLSchemaAttribute("comparison_type", "comparison_types", "The comparison type.  IE not equal, equal, etc.  IF value [comparison_type] cutoff_ or match_ we return TRUE.  Example (RMSDMetric) cutoff=4.0, comparison_type=lt.  We return true if the RMSD is less than 4.0. Required. Choieces are: \n" + utility::to_string( types));

	attlist + XMLSchemaAttribute("composite_action", xs_string, "If you have a composite metric this can be `any`, `all` or a specific composite value type.");

	attlist + XMLSchemaAttribute("cutoff", xsct_real, "Number to use to determine if filter passes or not for any RealMetric or CompositeRealMetric.  Set the comparison_type to indicate the behavior of this filter.");

	attlist + XMLSchemaAttribute("match", xs_string, "String to match on to determine if filter passes or not for any StringMetric or CompositeStringMetric.  Set the comparison type to indicate the behavior of this filter.");

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	SimpleMetricFactory::get_instance()->define_simple_metric_xml_schema( xsd );
	XMLSchemaSimpleSubelementList subelements;
	subelements.add_group_subelement( & SimpleMetricFactory::get_instance()->simple_metric_xml_schema_group_name );

	protocols::filters::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, class_name(), "Run a SimpleMetric (Real) as a filter.  Set the cutoff and comparison_type to control the behavior of this filter. ", attlist, subelements );
}

bool
SimpleMetricFilter::compare_metric(core::Real value) const {
	using namespace numeric;

	switch( type_ ){
	case( bogus ) :
		throw CREATE_EXCEPTION(utility::excn::BadInput, "A comparison_type must be set.  ");

	case ( eq ) :
		TR << value << " eq " << cutoff_ <<" ? "<< std::endl;
		if ( equal_by_epsilon(value,cutoff_,epsilon_) ) return true;
		break;

	case ( ne ) :
		TR << value << " ne " << cutoff_ <<" ? "<< std::endl;
		if ( ! (equal_by_epsilon(value,cutoff_,epsilon_) ) ) return true;
		break;

	case ( lt ) :
		TR << value << " lt " << cutoff_ <<" ? "<< std::endl;
		if ( value < cutoff_ ) return true;
		break;

	case ( gt ) :
		TR << value << " gt " << cutoff_ <<" ? "<< std::endl;
		if ( value > cutoff_ ) return true;
		break;

	case ( lt_or_eq ) :
		TR << value << " lt_or_eq " << cutoff_ <<" ? "<< std::endl;
		if ( value < cutoff_ ) return true;
		if ( equal_by_epsilon(value,cutoff_,epsilon_) ) return true;
		break;

	case ( gt_or_eq ) :
		TR << value << " gt_or_eq " << cutoff_ <<" ? "<< std::endl;
		if ( value > cutoff_ ) return true;
		if ( equal_by_epsilon(value,cutoff_,epsilon_) ) return true;
		break;

	}
	return false;
}

bool
SimpleMetricFilter::compare_metric(int value) const {
	using namespace numeric;

	int local_cutoff = int(cutoff_);
	switch( type_ ){
	case( bogus ) :
		throw CREATE_EXCEPTION(utility::excn::BadInput, "A comparison_type must be set.  ");

	case ( eq ) :
		TR << value << " eq " << local_cutoff <<" ? "<< std::endl;
		if ( value == local_cutoff ) return true;
		break;

	case ( ne ) :
		TR << value << " ne " << local_cutoff <<" ? "<< std::endl;
		if ( value != local_cutoff ) return true;
		break;

	case ( lt ) :
		TR << value << " lt " << local_cutoff <<" ? "<< std::endl;
		if ( value < local_cutoff ) return true;
		break;

	case ( gt ) :
		TR << value << " gt " << local_cutoff <<" ? "<< std::endl;
		if ( value > local_cutoff ) return true;
		break;

	case ( lt_or_eq ) :
		TR << value << " lt_or_eq " << local_cutoff <<" ? "<< std::endl;
		if ( value <= local_cutoff ) return true;
		break;

	case ( gt_or_eq ) :
		TR << value << " gt_or_eq " << local_cutoff <<" ? "<< std::endl;
		if ( value >= local_cutoff ) return true;
		break;

	}
	return false;
}

bool
SimpleMetricFilter::compare_metric( std::string const & value) const {

	if ( match_ == "" ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput, "For a StringMetric or CompositeStringMetric you must set the match string.");
	}

	TR << "Value: " << value << std::endl;
	switch( type_ ){
	case( bogus ) :
		utility_exit_with_message("A comparison_type must be set.  ");

	case ( eq ) :
		TR << value << " eq " << match_ <<" ? "<< std::endl;
		if ( value == match_ ) return true;
		break;

	case ( ne ) :
		TR << value << " ne " << match_ <<" ? "<< std::endl;
		if ( value != match_  ) return true;
		break;

	default :
		throw CREATE_EXCEPTION(utility::excn::BadInput, "SimpleMetricFilter: eq and ne are the only comparison_types that can be used with StringMetrics.");
	}
	return false;
}

template <class T>
bool SimpleMetricFilter::compare_composites( std::map<std::string, T > const & values) const {

	//std::cout << "Comparing Composites: "<<":"<<composite_action_<<":" << std::endl;
	if ( composite_action_ == "" ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput, "composite_action must be set in order to work with CompositeMetrics.");
	}

	if ( composite_action_ == "any" ) {
		for ( auto v_pair : values ) {
			if ( compare_metric( v_pair.second) ) {
				TR << v_pair.first << " passes set cutoff " << cutoff_ << std::endl;
				return true;
			}
		}
		return false;
	} else if ( composite_action_ == "all" ) {
		for ( auto v_pair : values ) {
			if ( !compare_metric( v_pair.second) ) {
				TR << v_pair.first << " did not pass cutuff."<< cutoff_ <<  "Filtering. " << std::endl;
				return false;
			}
		}
		return true;
	} else if ( values.count(composite_action_) ) {
		for ( auto v_pair : values ) {
			if ( v_pair.first == composite_action_ ) {
				bool pass = compare_metric( v_pair.second );
				TR << "metric type " << composite_action_ << " Pass? " << pass << std::endl;
				return pass;
			}
		}
	} else {
		TR << "Composite Action not understood :" << composite_action_ << ":" << std::endl;
		TR << "This can be `all`, `any`, or any of the following metric value types from the composite metric:" << std::endl;
		for ( auto metric_pair : values ) {
			TR << "  \t  " << metric_pair.first << " , ";
		}
		TR << std::endl;
		throw CREATE_EXCEPTION(utility::excn::BadInput, " Bad input for composite_action in SimpleMetricFilter. :"+composite_action_+":");
	}
	return false;
}

bool
SimpleMetricFilter::apply( core::pose::Pose const & pose) const
{
	using namespace numeric;

	if ( ! metric_ ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput, "SimpleMetricFilter requires a metric to be set.");
	}

	//This is just a triple check.  We should never ever be here.
	if ( metric_->type() == "RealMetric" ) {
		RealMetric const & r_metric = dynamic_cast<RealMetric const & >( *metric_ );
		core::Real value = r_metric.calculate(pose);
		bool pass = compare_metric( value );
		TR << "Filter passed: " << pass << std::endl;
		return pass;
	} else if ( metric_->type() == "StringMetric" ) {
		StringMetric const & r_metric = dynamic_cast<StringMetric const & >( *metric_ );
		std::string value = r_metric.calculate(pose);
		bool pass = compare_metric( value );
		TR << "Filter passed: " << pass << std::endl;
		return pass;
	} else if ( metric_->type() == "IntegerMetric" ) {
		IntegerMetric const & r_metric = dynamic_cast<IntegerMetric const & >( *metric_ );
		int value = r_metric.calculate(pose);
		bool pass = compare_metric( value );
		TR << "Filter passed: " << pass << std::endl;
		return pass;

	} else if ( metric_ ->type() == "CompositeRealMetric" ) {
		CompositeRealMetric const & r_metric = dynamic_cast<CompositeRealMetric const & >( *metric_ );
		std::map< std::string, core::Real > values = r_metric.calculate( pose );
		bool pass = compare_composites(values);
		TR << "Filter passed: " << pass << std::endl;
		return pass;
	} else if ( metric_ ->type() == "CompositeIntegerMetric" ) {
		CompositeIntegerMetric const & r_metric = dynamic_cast<CompositeIntegerMetric const & >( *metric_ );
		std::map< std::string, int > values = r_metric.calculate( pose );
		bool pass = compare_composites(values);
		TR << "Filter passed: " << pass << std::endl;
		return pass;
	} else if ( metric_->type() == "CompositeStringMetric" ) {
		CompositeStringMetric const & r_metric = dynamic_cast<CompositeStringMetric const & >( *metric_ );
		std::map< std::string, std::string > values = r_metric.calculate( pose );
		std::cout << "Comparing String Metrics " <<std::endl;
		bool pass = compare_composites(values);
		TR << "Filter passed: " << pass << std::endl;
		return pass;
	} else {
		utility_exit_with_message("SimpleMetric type not compatible. "+metric_->type());
	}
}

/////////////// Creator ///////////////

protocols::filters::FilterOP
SimpleMetricFilterCreator::create_filter() const
{
	return protocols::filters::FilterOP( new SimpleMetricFilter );
}

std::string
SimpleMetricFilterCreator::keyname() const
{
	return SimpleMetricFilter::class_name();
}

void SimpleMetricFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SimpleMetricFilter::provide_xml_schema( xsd );
}

} //protocols
} //simple_filters
