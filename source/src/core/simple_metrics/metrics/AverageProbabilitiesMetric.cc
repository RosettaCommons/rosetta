// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/AverageProbabilitiesMetric.cc
/// @brief A metric for averaging multiple PerResidueProbabilitiesMetrics
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)

// Unit headers
#include <core/simple_metrics/metrics/AverageProbabilitiesMetric.hh>
#include <core/simple_metrics/metrics/AverageProbabilitiesMetricCreator.hh>

// Core headers
#include <core/simple_metrics/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/pointer/memory.hh>

// std headers
#include <numeric>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "core.simple_metrics.metrics.AverageProbabilitiesMetric" );


namespace core {
namespace simple_metrics {
namespace metrics {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
AverageProbabilitiesMetric::AverageProbabilitiesMetric():
	core::simple_metrics::PerResidueProbabilitiesMetric()
{}

AverageProbabilitiesMetric::AverageProbabilitiesMetric( utility::vector1<core::simple_metrics::SimpleMetricCOP> metrics, utility::vector1< core::Real > weights ):
	core::simple_metrics::PerResidueProbabilitiesMetric()
{
	set_metric( metrics, weights );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
AverageProbabilitiesMetric::~AverageProbabilitiesMetric()= default;

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
AverageProbabilitiesMetric::AverageProbabilitiesMetric(AverageProbabilitiesMetric const &  ) = default;

core::simple_metrics::SimpleMetricOP
AverageProbabilitiesMetric::clone() const {
	return utility::pointer::make_shared< AverageProbabilitiesMetric >(*this );
}

std::string
AverageProbabilitiesMetric::name() const {
	return name_static();
}

std::string
AverageProbabilitiesMetric::name_static() {
	return "AverageProbabilitiesMetric";

}
std::string
AverageProbabilitiesMetric::metric() const {
	return name_static();
}

void
AverageProbabilitiesMetric::set_metric(utility::vector1<core::simple_metrics::SimpleMetricCOP> metrics, utility::vector1< core::Real > weights ){
	if ( metrics.size() != weights.size() ) {
		utility_exit_with_message("Number of metrics does not match number of weights!");
	}
	metrics_ = metrics;
	weights_ = weights;
}

void
AverageProbabilitiesMetric::set_use_cached_data(bool use_cache, std::string const & prefix, std::string const & suffix){
	use_cache_ = use_cache;
	cache_prefix_ = prefix;
	cache_suffix_ = suffix;
	if ( use_cache_ ) {
		TR << "Attempting to use cached data with set prefix/suffix:" <<prefix <<" "<<suffix << std::endl;
	}
}

void
AverageProbabilitiesMetric::set_fail_on_missing_cache(bool fail){
	fail_on_missing_cache_ = fail;
}

void
AverageProbabilitiesMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data )
{
	SimpleMetric::parse_base_tag( tag );

	std::string metrics_str = tag->getOption<std::string>("metrics");

	utility::vector1<std::string> metric_names = utility::string_split(metrics_str, ',');

	utility::vector1<core::simple_metrics::SimpleMetricCOP> metrics = core::simple_metrics::get_metrics_from_datamap_and_subtags(tag, data, "metrics" );

	std::string weights_str = tag->getOption<std::string>("weights", "");

	// if no weights are provided, set all to 1.0
	utility::vector1<core::Real> weights;
	if ( weights_str.empty() ) {
		core::Size num_metrics = metrics.size();
		for ( core::Size i = 0; i < num_metrics; ++i ) {
			weights.push_back( 1.0 );
		}
	} else {
		weights = utility::string_split<core::Real>(weights_str, ',', core::Real());
	}

	set_metric( metrics, weights );

	bool use_cache = tag->getOption<bool>("use_cached_data", false);
	std::string prefix;
	std::string suffix;

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
AverageProbabilitiesMetric::provide_xml_schema(utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;

	attlist + XMLSchemaAttribute::required_attribute( "metrics", xs_string, "A list of comma-seperated PerResidueProbabilitiesMetrics to calculate the average of." );

	attlist + XMLSchemaAttribute::attribute_w_default( "weights", xs_string, "A list of comma-seperated values to weight each metric with (Defaults to 1.0 for all). You need to provide as many weights as you provide metrics.", "");

	attlist + XMLSchemaAttribute::attribute_w_default( "use_cached_data",  xsct_rosetta_bool, "Use any data stored in the datacache that matches the set metrics name (and any prefix/suffix.)  Data is stored during a SimpleMetric's apply function, which is called during RunSimpleMetrics", "false");
	attlist + XMLSchemaAttribute("cache_prefix", xs_string, "Any prefix used during apply (RunSimpleMetrics), that we will match on if use_cache is true. Requires that all PerResidueProbabilitiesMetrics have the same prefix/suffix.");
	attlist + XMLSchemaAttribute("cache_suffix", xs_string, "Any suffix used during apply (RunSimpleMetrics), that we will match on if use_cache is true. Requires that all PerResidueProbabilitiesMetrics have the same prefix/suffix.");
	attlist + XMLSchemaAttribute::attribute_w_default("fail_on_missing_cache", xsct_rosetta_bool, "If use_cached_data is True and cache is not found, should we fail?", "true");

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"A metric for calculating the weighted average of multiple PerResidueProbabilitiesMetrics.", attlist);

}

std::map<core::Size, std::map<core::chemical::AA, core::Real>>
AverageProbabilitiesMetric::calculate(core::pose::Pose const &pose ) const {

	utility::vector1< std::map< core::Size, std::map< core::chemical::AA, core::Real >>> all_values;
	for ( auto const & metricCOP : metrics_ ) {

		if ( metricCOP->simple_metric_type() != "PerResidueProbabilitiesMetric" ) {
			utility_exit_with_message("AverageProbabilitiesMetric only works with PerResidueProbabilitiesMetrics!");
		}
		core::simple_metrics::PerResidueProbabilitiesMetricCOP res_metric = utility::pointer::dynamic_pointer_cast< core::simple_metrics::PerResidueProbabilitiesMetric const>( metricCOP );
		std::map< core::Size, std::map< core::chemical::AA, core::Real >> const values = res_metric->cached_calculate( pose, use_cache_, cache_prefix_, cache_suffix_, fail_on_missing_cache_ );
		all_values.push_back(values);
	}

	std::map<core::Size, std::map<core::chemical::AA, core::Real>> average_map = compute_average( all_values );

	return average_map;
}

std::map<core::Size, std::map<core::chemical::AA, core::Real>>
AverageProbabilitiesMetric::compute_average(utility::vector1<std::map<core::Size, std::map<core::chemical::AA, core::Real>>> const & all_values) const {

	std::map<core::Size, std::map<core::chemical::AA, core::Real>> averaged_values;

	// Calculate the total weight for normalization.
	core::Real total_weight = std::accumulate(weights_.begin(), weights_.end(), 0.0);

	// Iterate over each set of values together with its corresponding weight.
	for ( core::Size metric_idx = 1; metric_idx <= all_values.size(); ++metric_idx ) {
		auto const& values = all_values[metric_idx];
		core::Real weight = weights_[metric_idx];
		core::Real normalized_weight = weight / total_weight;
		TR << "WEIGHT IS: " << weight << " normalized weight is " << normalized_weight << std::endl;

		for ( auto const& pos_and_probs : values ) {
			core::Size const& position = pos_and_probs.first;
			for ( auto const& aa_and_prob : pos_and_probs.second ) {
				core::chemical::AA const& aa = aa_and_prob.first;
				core::Real const& prob = aa_and_prob.second;

				// If this amino acid hasn't been seen at this position, initialize its value to zero.
				if ( averaged_values[position].find(aa) == averaged_values[position].end() ) {
					averaged_values[position][aa] = 0.0;
				}

				// Update the averaged value by adding the weighted probability.
				averaged_values[position][aa] += prob * normalized_weight;
			}
		}
	}
	return averaged_values;
}

void
AverageProbabilitiesMetricCreator::provide_xml_schema(utility::tag::XMLSchemaDefinition & xsd ) const {
	AverageProbabilitiesMetric::provide_xml_schema(xsd );
}

std::string
AverageProbabilitiesMetricCreator::keyname() const {
	return AverageProbabilitiesMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
AverageProbabilitiesMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< AverageProbabilitiesMetric >();
}

void
AverageProbabilitiesMetric::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"AverageProbabilitiesMetric", basic::citation_manager::CitedModuleType::SimpleMetric,
		"Moritz Ertelt",
		"University of Leipzig",
		"moritz.ertelt@gmail.com",
		""
		)
	);
}

} //metrics
} //simple_metrics
} //core


#ifdef    SERIALIZATION

template< class Archive >
void
core::simple_metrics::metrics::AverageProbabilitiesMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::PerResidueProbabilitiesMetric>( this ) );
	arc( CEREAL_NVP( metrics_ ) );
	arc( CEREAL_NVP( use_cache_ ) );
	arc( CEREAL_NVP( cache_prefix_ ) );
	arc( CEREAL_NVP( cache_suffix_ ) );
	arc( CEREAL_NVP( fail_on_missing_cache_ ) );
	arc( CEREAL_NVP ( weights_ ) );

}

template< class Archive >
void
core::simple_metrics::metrics::AverageProbabilitiesMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::PerResidueProbabilitiesMetric >( this ) );

	utility::vector1<core::simple_metrics::SimpleMetricCOP> local_metrics;
	arc( local_metrics); // PerResidueProbabilitiesMetricCOP
	metrics_ = local_metrics; // copy the non-const pointer(s) into the const pointer(s)

	arc( metrics_ );
	arc( use_cache_ );
	arc( cache_prefix_ );
	arc( cache_suffix_ );
	arc( fail_on_missing_cache_ );
	arc( weights_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::simple_metrics::metrics::AverageProbabilitiesMetric );
CEREAL_REGISTER_TYPE( core::simple_metrics::metrics::AverageProbabilitiesMetric )

CEREAL_REGISTER_DYNAMIC_INIT( core_simple_metrics_metrics_AverageProbabilitiesMetric )
#endif // SERIALIZATION





