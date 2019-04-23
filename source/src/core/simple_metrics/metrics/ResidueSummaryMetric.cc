// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/ResidueSummaryMetric.cc
/// @brief A SimpleMetric that takes a PerResidueMetric and calculates different summaries of the overall data.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <core/simple_metrics/metrics/ResidueSummaryMetric.hh>
#include <core/simple_metrics/PerResidueRealMetric.hh>
#include <core/simple_metrics/simple_metric_creators.hh>

// Core headers
#include <core/simple_metrics/RealMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/util.hh>
#include <utility/string_util.hh>
#include <numeric/util.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "core.simple_metrics.metrics.ResidueSummaryMetric" );


namespace core {
namespace simple_metrics {
namespace metrics {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
ResidueSummaryMetric::ResidueSummaryMetric():
	core::simple_metrics::RealMetric()
{}

ResidueSummaryMetric::ResidueSummaryMetric(PerResidueRealMetricCOP metric):
	core::simple_metrics::RealMetric()
{
	set_metric(metric);
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
ResidueSummaryMetric::~ResidueSummaryMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
ResidueSummaryMetric::ResidueSummaryMetric( ResidueSummaryMetric const &  ) = default;

core::simple_metrics::SimpleMetricOP
ResidueSummaryMetric::clone() const {
	return core::simple_metrics::SimpleMetricOP(new ResidueSummaryMetric( *this ) );

}

std::string
ResidueSummaryMetric::name() const {
	return name_static();
}

std::string
ResidueSummaryMetric::name_static() {
	return "ResidueSummaryMetric";

}
std::string
ResidueSummaryMetric::metric() const {
	return "res_summary";
}

void
ResidueSummaryMetric::set_action_value(core::Real cutoff){
	cutoff_ = cutoff;
}

void
ResidueSummaryMetric::set_epsilon( core::Real epsilon ){
	epsilon_ = epsilon;
}

void
ResidueSummaryMetric::set_action(core::simple_metrics::metrics::summary_type sum_type){
	action_ = sum_type;
}

void
ResidueSummaryMetric::set_metric(core::simple_metrics::PerResidueRealMetricCOP metric){
	metric_ = metric;
}

void
ResidueSummaryMetric::set_use_cached_data(bool use_cache, std::string prefix, std::string suffix){
	use_cache_ = use_cache;
	cache_prefix_ = prefix;
	cache_suffix_ = suffix;
	if ( use_cache_ ) {
		TR << "Attempting to use cached data with set prefix/suffix:" <<prefix <<" "<<suffix << std::endl;
	}
}

void
ResidueSummaryMetric::set_fail_on_missing_cache(bool fail){
	fail_on_missing_cache_ = fail;
}

void
ResidueSummaryMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data )
{
	using namespace utility::pointer;

	SimpleMetric::parse_base_tag( tag );

	SimpleMetricCOP metric = get_metric_from_datamap_and_subtags(tag, data, "metric");

	if ( metric->simple_metric_type() != "PerResidueRealMetric" ) {
		utility_exit_with_message("ResidueSummaryMetric only works with PerResidueRealMetrics!");
	}

	PerResidueRealMetricCOP res_metric = dynamic_pointer_cast< PerResidueRealMetric const>( metric );

	set_metric(res_metric);

	set_action( summary_string_to_type.at(tag->getOption< std::string >("action")));
	set_epsilon( tag->getOption< core::Real >("epsilon", epsilon_));
	set_action_value( tag->getOption< core::Real >("action_value", 0));


	bool use_cache = tag->getOption<bool>("use_cached_data", false);
	std::string prefix= "";
	std::string suffix="";

	if ( tag->hasOption("cache_prefix") ) {
		prefix = tag->getOption< std::string >("cache_prefix");
	}
	if ( tag->hasOption("cache_suffix") ) {
		suffix = tag->getOption< std::string >("cache_suffix");
	}

	set_use_cached_data(use_cache, prefix, suffix);
	set_fail_on_missing_cache(tag->getOption< bool>("fail_on_missing_cache", fail_on_missing_cache_));

}

void
ResidueSummaryMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;

	std::string epsilon_description =
		"Significant digits used for equal operators\n"
		"\n"
		" IE:\n"
		"  The absolute value of the difference between two numbers, below which they will be considered to be equal.\n"
		"  Used when comparing numbers if action = n_res_eq or n_res_ne.";

	attlist + XMLSchemaAttribute::attribute_w_default("epsilon", xsct_real , epsilon_description, ".0001");
	attlist + XMLSchemaAttribute( "metric", xs_string, "The PerResidueRealMetric that we will summarize ");


	utility::vector1< std::string > types = get_summary_type_strings();
	utility::tag::add_schema_restrictions_for_strings( xsd, "summary_types", types);

	std::string const action_docs = "Summary type that this metric should calculate.  Current choices are:\n" + utility::to_string(types);

	attlist + XMLSchemaAttribute::attribute_w_default("action", "summary_types", action_docs, "mean");

	attlist + XMLSchemaAttribute("action_value", xsct_real, "Number to use for any summary metrics that calculate N residues based on a value. ");

	//attributes_for_parse_residue_selector( attlist, "residue_selector",
	// "Selector specifying residues." );

	attlist + XMLSchemaAttribute::attribute_w_default( "use_cached_data",  xsct_rosetta_bool, "Use any data stored in the datacache that matches the set metrics name (and any prefix/suffix.)  Data is stored during a SimpleMetric's apply function, which is called during RunSimpleMetrics", "false");
	attlist + XMLSchemaAttribute("cache_prefix", xs_string, "Any prefix used during apply (RunSimpleMetrics), that we will match on if use_cache is true");
	attlist + XMLSchemaAttribute("cache_suffix", xs_string, "Any suffix used during apply (RunSimpleMetrics), that we will match on if use_cache is true");
	attlist + XMLSchemaAttribute::attribute_w_default("fail_on_missing_cache", xsct_rosetta_bool, "If use_cached_data is True and cache is not found, should we fail?", "true");

	std::string description =
		"Author: Jared Adolf-Bryfogle (jadolfbr@gmail.com)\n"
		"A SimpleMetric that takes a PerResidueRealMetric and calculates different summaries of the overall data.\n"
		"  This metric can calculate means, totals, or the number of residues (n_res) matching certain criteria. \n"
		"  Useful for summarizing metrics or using more complex functionality in the SimpleMetricFilter\n"
		"\n"
		"  Be sure to set a custom_type to label the summary type in which you are calculating!.";

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		description, attlist);
}

core::Real
ResidueSummaryMetric::calculate(const core::pose::Pose & pose) const {

	if ( metric_ == nullptr ) {
		utility_exit_with_message("ResidueSummaryMetric: This metric requires a PerResidueRealMetric to calculate its summary");
	}


	std::map< core::Size, core::Real > const values = metric_->cached_calculate( pose, use_cache_, cache_prefix_, cache_suffix_, fail_on_missing_cache_ );

	core::Size n_res = 0;
	core::Real sum_of_values = 0.0;

	switch( action_ ){
	case( bogus ) :
		utility_exit_with_message("A summary_type must be set.  ");

	case ( sum ) :
		TR << "Returning the sum of the data" << std::endl;

		sum_of_values = 0.0;
		for ( auto const & res_value: values ) {
			sum_of_values+= res_value.second;
		}
		return sum_of_values;

	case ( mean ) :
		TR << "Returning the mean of the data" << std::endl;

		sum_of_values = 0.0;
		for ( auto const & res_value: values ) {
			sum_of_values += res_value.second;
		}
		return sum_of_values/values.size();

	case ( n_res_eq ) :
		TR << "Returning n res equal to "<< cutoff_ << std::endl;

		n_res = 0;
		for ( auto const & res_value: values ) {
			if ( numeric::equal_by_epsilon( res_value.second,cutoff_,epsilon_ ) ) {
				n_res += 1;
			}
		}
		return n_res;

	case ( n_res_ne) :
		TR << "Returning n res not equal to " << cutoff_ << std::endl;

		n_res = 0;
		for ( auto const & res_value: values ) {
			if ( ! numeric::equal_by_epsilon( res_value.second,cutoff_,epsilon_ ) ) {
				n_res += 1;
			}
		}

		return n_res;

	case ( n_res_lt ) :
		TR << "Returning n res less than "<< cutoff_ << std::endl;

		n_res = 0;
		for ( auto const & res_value : values ) {
			if ( res_value.second < cutoff_ ) {
				n_res += 1;
			}
		}
		return n_res;
		break;

	case ( n_res_lt_or_eq ) :
		TR << "Returning n res less than or equal to "<< cutoff_ << std::endl;

		n_res = 0;
		for ( auto const & res_value : values ) {
			if ( res_value.second <= cutoff_ ) {
				n_res += 1;
			}
		}
		return n_res;

	case ( n_res_gt ) :
		TR << "Returning n res greater than "<< cutoff_ << std::endl;

		n_res = 0;
		for ( auto const & res_value : values ) {
			if ( res_value.second > cutoff_ ) {
				n_res += 1;
			}
		}
		return n_res;

	case ( n_res_gt_or_eq) :
		TR << "Returning n res greater than or equal to "<< cutoff_ << std::endl;

		n_res = 0;
		for ( auto const & res_value : values ) {
			if ( res_value.second >= cutoff_ ) {
				n_res += 1;
			}
		}
		return n_res;

	}

	return 0;
}



void
ResidueSummaryMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ResidueSummaryMetric::provide_xml_schema( xsd );
}

std::string
ResidueSummaryMetricCreator::keyname() const {
	return ResidueSummaryMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
ResidueSummaryMetricCreator::create_simple_metric() const {
	return core::simple_metrics::SimpleMetricOP( new ResidueSummaryMetric );

}

utility::vector1< std::string >
get_summary_type_strings(){
	utility::vector1< std::string > type_list;
	for ( auto type_pair : summary_string_to_type ) {
		type_list.push_back(type_pair.first);
	}
	return type_list;
}





} //core
} //simple_metrics
} //metrics


#ifdef    SERIALIZATION



template< class Archive >
void
core::simple_metrics::metrics::ResidueSummaryMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::RealMetric>( this ) );
	arc( CEREAL_NVP( metric_ ) );
	arc( CEREAL_NVP(cutoff_) );
	arc( CEREAL_NVP( epsilon_ ) );
	arc( CEREAL_NVP( action_ ) );
	arc( CEREAL_NVP( use_cache_ ) );
	arc( CEREAL_NVP( cache_prefix_ ) );
	arc( CEREAL_NVP( cache_suffix_ ) );
	arc( CEREAL_NVP( fail_on_missing_cache_ ) );
}

template< class Archive >
void
core::simple_metrics::metrics::ResidueSummaryMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::RealMetric >( this ) );

	std::shared_ptr< core::simple_metrics::PerResidueRealMetric > local_metric;
	arc( local_metric); // PerResidueRealMetricCOP
	metric_ = local_metric; // copy the non-const pointer(s) into the const pointer(s)

	arc( metric_ );
	arc( cutoff_ );
	arc( epsilon_ );
	arc( action_ );
	arc( use_cache_ );
	arc( cache_prefix_ );
	arc( cache_suffix_ );
	arc( fail_on_missing_cache_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::simple_metrics::metrics::ResidueSummaryMetric );
CEREAL_REGISTER_TYPE( core::simple_metrics::metrics::ResidueSummaryMetric )

CEREAL_REGISTER_DYNAMIC_INIT( core_simple_metrics_metrics_ResidueSummaryMetric )
#endif // SERIALIZATION




