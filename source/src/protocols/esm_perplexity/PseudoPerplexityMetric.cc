// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/esm_perplexity/PseudoPerplexityMetric.cc
/// @brief A class for calculating the pseudo-perplexity from a given PerResidueProbabilitiesMetric.
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)
/// @note This has been adopted from how the ResidueSummaryMetric from Jared works.

// Unit headers
#include <protocols/esm_perplexity/PseudoPerplexityMetric.hh>
#include <protocols/esm_perplexity/PseudoPerplexityMetricCreator.hh>
#include <protocols/esm_perplexity/EsmPerplexityTensorflowProtocol.hh>

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

static basic::Tracer TR( "protocols.esm_perplexity.PseudoPerplexityMetric" );


namespace protocols {
namespace esm_perplexity {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
PseudoPerplexityMetric::PseudoPerplexityMetric():
	core::simple_metrics::RealMetric()
{}

PseudoPerplexityMetric::PseudoPerplexityMetric( core::simple_metrics::PerResidueProbabilitiesMetricCOP metric ):
	core::simple_metrics::RealMetric()
{
	set_metric( metric );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
PseudoPerplexityMetric::~PseudoPerplexityMetric()= default;

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
PseudoPerplexityMetric::PseudoPerplexityMetric(PseudoPerplexityMetric const &  ) = default;

core::simple_metrics::SimpleMetricOP
PseudoPerplexityMetric::clone() const {
	return utility::pointer::make_shared< PseudoPerplexityMetric >(*this );
}

std::string
PseudoPerplexityMetric::name() const {
	return name_static();
}

std::string
PseudoPerplexityMetric::name_static() {
	return "PseudoPerplexityMetric";

}
std::string
PseudoPerplexityMetric::metric() const {
	return name_static();
}

void
PseudoPerplexityMetric::set_metric(core::simple_metrics::PerResidueProbabilitiesMetricCOP metric){
	metric_ = metric;
}

void
PseudoPerplexityMetric::set_use_cached_data(bool use_cache, std::string const & prefix, std::string const & suffix){
	use_cache_ = use_cache;
	cache_prefix_ = prefix;
	cache_suffix_ = suffix;
	if ( use_cache_ ) {
		TR << "Attempting to use cached data with set prefix/suffix:" <<prefix <<" "<<suffix << std::endl;
	}
}

void
PseudoPerplexityMetric::set_fail_on_missing_cache(bool fail){
	fail_on_missing_cache_ = fail;
}

/// @brief Function to return the (pseudo-)perplexity from a map of probabilities
core::Real
PseudoPerplexityMetric::compute_perplexity(
	core::pose::Pose const & pose,
	std::map< core::Size, std::map< core::chemical::AA, core::Real>> const & values
) {
	core::Real log_probabilities_sum = 0;
	for ( auto const & position_pair : values ) {
		core::chemical::AA aa_type = pose.residue( position_pair.first ).aa();
		if ( aa_type == core::chemical::aa_unk ) {
			utility_exit_with_message("Residue type is unknown, this might lead to ambiguity and shouldn't have happened in the first place.");
		}
		core::Real aa_probability = position_pair.second.at( aa_type );
		// if the probability is zero add a small constant to it to prevent -inf
		if ( aa_probability == 0 || aa_probability < 0.00001 ) {
			aa_probability += 0.00001;
		}
		if ( std::isnan(aa_probability) || std::isinf(aa_probability) ) {
			aa_probability = 0.00001;
		}
		log_probabilities_sum += std::log( aa_probability );
	}
	core::Real perplexity = std::exp( -log_probabilities_sum/static_cast<core::Real>(values.size()));
	return perplexity;
}

void
PseudoPerplexityMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data )
{
	SimpleMetric::parse_base_tag( tag );

	core::simple_metrics::SimpleMetricCOP metric = core::simple_metrics::get_metric_from_datamap_and_subtags(tag, data, "metric");

	if ( metric->simple_metric_type() != "PerResidueProbabilitiesMetric" ) {
		utility_exit_with_message("PseudoPerplexityMetric only works with PerResidueProbabilitiesMetrics!");
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
PseudoPerplexityMetric::provide_xml_schema(utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;

	attlist + XMLSchemaAttribute::required_attribute( "metric", xs_string, "A PerResidueProbabilitiesMetric to calculate the pseudo-perplexity from." );

	attlist + XMLSchemaAttribute::attribute_w_default( "use_cached_data",  xsct_rosetta_bool, "Use any data stored in the datacache that matches the set metrics name (and any prefix/suffix.)  Data is stored during a SimpleMetric's apply function, which is called during RunSimpleMetrics", "false");
	attlist + XMLSchemaAttribute("cache_prefix", xs_string, "Any prefix used during apply (RunSimpleMetrics), that we will match on if use_cache is true");
	attlist + XMLSchemaAttribute("cache_suffix", xs_string, "Any suffix used during apply (RunSimpleMetrics), that we will match on if use_cache is true");
	attlist + XMLSchemaAttribute::attribute_w_default("fail_on_missing_cache", xsct_rosetta_bool, "If use_cached_data is True and cache is not found, should we fail?", "true");

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"A metric for estimating the likeliness of a sequence given some predicted probabilities.", attlist);

}

core::Real
PseudoPerplexityMetric::calculate(core::pose::Pose const &pose ) const {

	if ( metric_ == nullptr ) {
		utility_exit_with_message("PseudoPerplexityMetric: This metric requires a PerResidueProbabilitiesMetric!");
	}

	std::map< core::Size, std::map< core::chemical::AA, core::Real >> const values = metric_->cached_calculate( pose, use_cache_, cache_prefix_, cache_suffix_, fail_on_missing_cache_ );

	core::Real pseudo_perplexity = compute_perplexity( pose, values );
	return pseudo_perplexity;
}

void
PseudoPerplexityMetricCreator::provide_xml_schema(utility::tag::XMLSchemaDefinition & xsd ) const {
	PseudoPerplexityMetric::provide_xml_schema(xsd );
}

std::string
PseudoPerplexityMetricCreator::keyname() const {
	return PseudoPerplexityMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
PseudoPerplexityMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< PseudoPerplexityMetric >();
}

void
PseudoPerplexityMetric::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"PseudoPerplexityMetric", basic::citation_manager::CitedModuleType::SimpleMetric,
		"Moritz Ertelt",
		"University of Leipzig",
		"moritz.ertelt@gmail.com",
		""
		)
	);
}

} //esm_perplexity
} //protocols


#ifdef    SERIALIZATION

template< class Archive >
void
protocols::esm_perplexity::PseudoPerplexityMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::RealMetric>( this ) );
	arc( CEREAL_NVP( metric_ ) );
	arc( CEREAL_NVP( use_cache_ ) );
	arc( CEREAL_NVP( cache_prefix_ ) );
	arc( CEREAL_NVP( cache_suffix_ ) );
	arc( CEREAL_NVP( fail_on_missing_cache_ ) );

}

template< class Archive >
void
protocols::esm_perplexity::PseudoPerplexityMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::RealMetric >( this ) );

	std::shared_ptr< core::simple_metrics::PerResidueProbabilitiesMetric > local_metric;
	arc( local_metric); // PerResidueRealMetricCOP
	metric_ = local_metric; // copy the non-const pointer(s) into the const pointer(s)

	arc( metric_ );
	arc( use_cache_ );
	arc( cache_prefix_ );
	arc( cache_suffix_ );
	arc( fail_on_missing_cache_ );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::esm_perplexity::PseudoPerplexityMetric );
CEREAL_REGISTER_TYPE( protocols::esm_perplexity::PseudoPerplexityMetric )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_esm_perplexity_PseudoPerplexityMetric )
#endif // SERIALIZATION





