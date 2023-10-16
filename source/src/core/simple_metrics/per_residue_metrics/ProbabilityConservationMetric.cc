// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/ProbabilityConservationMetric.cc
/// @brief A class for calculating the conservation of a position given some predicted probabilities (using the relative Shannon Entropy).
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)

// Unit headers
#include <core/simple_metrics/per_residue_metrics/ProbabilityConservationMetric.hh>
#include <core/simple_metrics/per_residue_metrics/ProbabilityConservationMetricCreator.hh>

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

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "core.simple_metrics.per_residue_metrics.ProbabilityConservationMetric" );


namespace core {
namespace simple_metrics {
namespace per_residue_metrics {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
ProbabilityConservationMetric::ProbabilityConservationMetric():
	core::simple_metrics::PerResidueRealMetric()
{}

ProbabilityConservationMetric::ProbabilityConservationMetric( core::simple_metrics::PerResidueProbabilitiesMetricCOP metric ):
	core::simple_metrics::PerResidueRealMetric()
{
	set_metric( metric );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
ProbabilityConservationMetric::~ProbabilityConservationMetric()= default;

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
ProbabilityConservationMetric::ProbabilityConservationMetric(ProbabilityConservationMetric const &  ) = default;

core::simple_metrics::SimpleMetricOP
ProbabilityConservationMetric::clone() const {
	return utility::pointer::make_shared< ProbabilityConservationMetric >(*this );
}

std::string
ProbabilityConservationMetric::name() const {
	return name_static();
}

std::string
ProbabilityConservationMetric::name_static() {
	return "ProbabilityConservationMetric";

}
std::string
ProbabilityConservationMetric::metric() const {
	return name_static();
}

void
ProbabilityConservationMetric::set_metric(core::simple_metrics::PerResidueProbabilitiesMetricCOP metric){
	metric_ = metric;
}

void
ProbabilityConservationMetric::set_use_cached_data(bool use_cache, std::string const & prefix, std::string const & suffix){
	use_cache_ = use_cache;
	cache_prefix_ = prefix;
	cache_suffix_ = suffix;
	if ( use_cache_ ) {
		TR << "Attempting to use cached data with set prefix/suffix:" <<prefix <<" "<<suffix << std::endl;
	}
}

void
ProbabilityConservationMetric::set_fail_on_missing_cache(bool fail){
	fail_on_missing_cache_ = fail;
}

std::map<core::Size, core::Real>
ProbabilityConservationMetric::compute_conservation( std::map<core::Size, std::map<core::chemical::AA, core::Real>> const & values) const {

	std::map<core::Size, core::Real> conservation_values;
	core::Real max_entropy;

	for ( auto const & position_and_probs : values ) {
		core::Size const & position = position_and_probs.first;
		const std::map<core::chemical::AA, core::Real>& aa_probs = position_and_probs.second;
		max_entropy = std::log( aa_probs.size() ); // assuming equal likeliness of amino acids

		core::Real entropy = 0.0;

		for ( auto const & aa_and_prob : aa_probs ) {
			core::Real const & prob = aa_and_prob.second;
			if ( prob > 0 ) {
				entropy -= prob * std::log(prob);
			}
		}

		core::Real relative_entropy = entropy / max_entropy;
		conservation_values[position] = 1.0 - relative_entropy;  // Convert to conservation score
	}

	return conservation_values;
}

void
ProbabilityConservationMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data )
{
	SimpleMetric::parse_base_tag( tag );

	core::simple_metrics::SimpleMetricCOP metric = core::simple_metrics::get_metric_from_datamap_and_subtags(tag, data, "metric");

	if ( metric->simple_metric_type() != "PerResidueProbabilitiesMetric" ) {
		utility_exit_with_message("ProbabilityConservationMetric only works with PerResidueProbabilitiesMetrics!");
	}

	core::simple_metrics::PerResidueProbabilitiesMetricCOP res_metric = utility::pointer::dynamic_pointer_cast< core::simple_metrics::PerResidueProbabilitiesMetric const>( metric );

	set_metric( res_metric );

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
ProbabilityConservationMetric::provide_xml_schema(utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;

	attlist + XMLSchemaAttribute::required_attribute( "metric", xs_string, "A PerResidueProbabilitiesMetric to calculate the conservation of residues." );

	attlist + XMLSchemaAttribute::attribute_w_default( "use_cached_data",  xsct_rosetta_bool, "Use any data stored in the datacache that matches the set metrics name (and any prefix/suffix.)  Data is stored during a SimpleMetric's apply function, which is called during RunSimpleMetrics", "false");
	attlist + XMLSchemaAttribute("cache_prefix", xs_string, "Any prefix used during apply (RunSimpleMetrics), that we will match on if use_cache is true");
	attlist + XMLSchemaAttribute("cache_suffix", xs_string, "Any suffix used during apply (RunSimpleMetrics), that we will match on if use_cache is true");
	attlist + XMLSchemaAttribute::attribute_w_default("fail_on_missing_cache", xsct_rosetta_bool, "If use_cached_data is True and cache is not found, should we fail?", "true");

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"A PerResidueRealMetric for calculating the conservation of a position given some predicted probabilities (using the relative Shannon Entropy).", attlist);

}

std::map< core::Size, core::Real >
ProbabilityConservationMetric::calculate(core::pose::Pose const &pose ) const {

	if ( metric_ == nullptr ) {
		utility_exit_with_message("ProbabilityConservationMetric: This metric requires a PerResidueProbabilitiesMetric!");
	}

	std::map< core::Size, std::map< core::chemical::AA, core::Real >> const values = metric_->cached_calculate( pose, use_cache_, cache_prefix_, cache_suffix_, fail_on_missing_cache_ );

	std::map< core::Size, core::Real > conservation_map = compute_conservation( values );

	return conservation_map;
}



void
ProbabilityConservationMetricCreator::provide_xml_schema(utility::tag::XMLSchemaDefinition & xsd ) const {
	ProbabilityConservationMetric::provide_xml_schema(xsd );
}

std::string
ProbabilityConservationMetricCreator::keyname() const {
	return ProbabilityConservationMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
ProbabilityConservationMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< ProbabilityConservationMetric >();
}

void
ProbabilityConservationMetric::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"ProbabilityConservationMetric", basic::citation_manager::CitedModuleType::SimpleMetric,
		"Moritz Ertelt",
		"University of Leipzig",
		"moritz.ertelt@gmail.com",
		""
		)
	);
}

} //per_residue_metrics
} //simple_metrics
} //core


#ifdef    SERIALIZATION

template< class Archive >
void
core::simple_metrics::per_residue_metrics::ProbabilityConservationMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::PerResidueRealMetric>( this ) );
	arc( CEREAL_NVP( metric_ ) );
	arc( CEREAL_NVP( use_cache_ ) );
	arc( CEREAL_NVP( cache_prefix_ ) );
	arc( CEREAL_NVP( cache_suffix_ ) );
	arc( CEREAL_NVP( fail_on_missing_cache_ ) );

}

template< class Archive >
void
core::simple_metrics::per_residue_metrics::ProbabilityConservationMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::PerResidueRealMetric >( this ) );

	std::shared_ptr< core::simple_metrics::PerResidueProbabilitiesMetric > local_metric;
	arc( local_metric); // PerResidueRealMetricCOP
	metric_ = local_metric; // copy the non-const pointer(s) into the const pointer(s)

	arc( metric_ );
	arc( use_cache_ );
	arc( cache_prefix_ );
	arc( cache_suffix_ );
	arc( fail_on_missing_cache_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::simple_metrics::per_residue_metrics::ProbabilityConservationMetric );
CEREAL_REGISTER_TYPE( core::simple_metrics::per_residue_metrics::ProbabilityConservationMetric )

CEREAL_REGISTER_DYNAMIC_INIT( core_simple_metrics_per_residue_metrics_ProbabilityConservationMetric )
#endif // SERIALIZATION





