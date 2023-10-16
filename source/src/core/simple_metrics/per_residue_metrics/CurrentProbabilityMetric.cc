// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/CurrentProbabilityMetric.cc
/// @brief A class for for returning just the probability of the amino acid currently present in the pose from a PerResidueProbabilitiesMetric.
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)

// Unit headers
#include <core/simple_metrics/per_residue_metrics/CurrentProbabilityMetric.hh>
#include <core/simple_metrics/per_residue_metrics/CurrentProbabilityMetricCreator.hh>

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

static basic::Tracer TR( "core.simple_metrics.per_residue_metrics.CurrentProbabilityMetric" );

namespace core {
namespace simple_metrics {
namespace per_residue_metrics {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
CurrentProbabilityMetric::CurrentProbabilityMetric():
	core::simple_metrics::PerResidueRealMetric()
{}

CurrentProbabilityMetric::CurrentProbabilityMetric(core::simple_metrics::PerResidueProbabilitiesMetricCOP metric ):
	core::simple_metrics::PerResidueRealMetric()
{
	set_metric( metric );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
CurrentProbabilityMetric::~CurrentProbabilityMetric()= default;

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
CurrentProbabilityMetric::CurrentProbabilityMetric(CurrentProbabilityMetric const &  ) = default;

core::simple_metrics::SimpleMetricOP
CurrentProbabilityMetric::clone() const {
	return utility::pointer::make_shared< CurrentProbabilityMetric >(*this );
}

std::string
CurrentProbabilityMetric::name() const {
	return name_static();
}

std::string
CurrentProbabilityMetric::name_static() {
	return "CurrentProbabilityMetric";

}
std::string
CurrentProbabilityMetric::metric() const {
	return name_static();
}

void
CurrentProbabilityMetric::set_metric(core::simple_metrics::PerResidueProbabilitiesMetricCOP metric){
	metric_ = metric;
}

void
CurrentProbabilityMetric::set_use_cached_data(bool use_cache, std::string const & prefix, std::string const & suffix){
	use_cache_ = use_cache;
	cache_prefix_ = prefix;
	cache_suffix_ = suffix;
	if ( use_cache_ ) {
		TR << "Attempting to use cached data with set prefix/suffix:" <<prefix <<" "<<suffix << std::endl;
	}
}

void
CurrentProbabilityMetric::set_fail_on_missing_cache(bool fail){
	fail_on_missing_cache_ = fail;
}

std::map<core::Size, core::Real>
CurrentProbabilityMetric::get_current_probabilities(std::map<core::Size, std::map<core::chemical::AA, core::Real> > const & values,
	core::pose::Pose const & pose) {

	std::map<core::Size, core::Real> current_probs;

	for ( auto const& position_and_probs : values ) {
		core::Size const& position = position_and_probs.first;
		std::map<core::chemical::AA, core::Real> const& aa_probs = position_and_probs.second;

		// Get the current AA from the pose
		core::chemical::AA current_aa = pose.residue(position).aa();

		// Get its probability
		auto it = aa_probs.find(current_aa);
		if ( it != aa_probs.end() ) {  // Make sure the amino acid is in the map
			core::Real current_prob = it->second;
			current_probs[position] = current_prob;
		}
	}
	return current_probs;
}

void
CurrentProbabilityMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data )
{
	SimpleMetric::parse_base_tag( tag );

	core::simple_metrics::SimpleMetricCOP metric = core::simple_metrics::get_metric_from_datamap_and_subtags(tag, data, "metric");

	if ( metric->simple_metric_type() != "PerResidueProbabilitiesMetric" ) {
		utility_exit_with_message("CurrentProbabilityMetric only works with PerResidueProbabilitiesMetrics!");
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
CurrentProbabilityMetric::provide_xml_schema(utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;

	attlist + XMLSchemaAttribute::required_attribute( "metric", xs_string, "A PerResidueProbabilitiesMetric to lookup the probability of the current amino acid in the pose from." );

	attlist + XMLSchemaAttribute::attribute_w_default( "use_cached_data",  xsct_rosetta_bool, "Use any data stored in the datacache that matches the set metrics name (and any prefix/suffix.)  Data is stored during a SimpleMetric's apply function, which is called during RunSimpleMetrics", "false");
	attlist + XMLSchemaAttribute("cache_prefix", xs_string, "Any prefix used during apply (RunSimpleMetrics), that we will match on if use_cache is true");
	attlist + XMLSchemaAttribute("cache_suffix", xs_string, "Any suffix used during apply (RunSimpleMetrics), that we will match on if use_cache is true");
	attlist + XMLSchemaAttribute::attribute_w_default("fail_on_missing_cache", xsct_rosetta_bool, "If use_cached_data is True and cache is not found, should we fail?", "true");

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"A PerResidueRealMetric for returning just the probability of the amino acid currently present in the pose from a PerResidueProbabilitiesMetric.", attlist);

}

std::map< core::Size, core::Real >
CurrentProbabilityMetric::calculate(core::pose::Pose const &pose ) const {

	if ( metric_ == nullptr ) {
		utility_exit_with_message("CurrentProbabilityMetric: This metric requires a PerResidueProbabilitiesMetric!");
	}

	std::map< core::Size, std::map< core::chemical::AA, core::Real >> const values = metric_->cached_calculate( pose, use_cache_, cache_prefix_, cache_suffix_, fail_on_missing_cache_ );

	std::map< core::Size, core::Real > probability_map = get_current_probabilities( values, pose );

	return probability_map;
}

void
CurrentProbabilityMetricCreator::provide_xml_schema(utility::tag::XMLSchemaDefinition & xsd ) const {
	CurrentProbabilityMetric::provide_xml_schema(xsd );
}

std::string
CurrentProbabilityMetricCreator::keyname() const {
	return CurrentProbabilityMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
CurrentProbabilityMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< CurrentProbabilityMetric >();
}

void
CurrentProbabilityMetric::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"CurrentProbabilityMetric", basic::citation_manager::CitedModuleType::SimpleMetric,
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
core::simple_metrics::per_residue_metrics::CurrentProbabilityMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::PerResidueRealMetric>( this ) );
	arc( CEREAL_NVP( metric_ ) );
	arc( CEREAL_NVP( use_cache_ ) );
	arc( CEREAL_NVP( cache_prefix_ ) );
	arc( CEREAL_NVP( cache_suffix_ ) );
	arc( CEREAL_NVP( fail_on_missing_cache_ ) );

}

template< class Archive >
void
core::simple_metrics::per_residue_metrics::CurrentProbabilityMetric::load( Archive & arc ) {
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

SAVE_AND_LOAD_SERIALIZABLE( core::simple_metrics::per_residue_metrics::CurrentProbabilityMetric );
CEREAL_REGISTER_TYPE( core::simple_metrics::per_residue_metrics::CurrentProbabilityMetric )

CEREAL_REGISTER_DYNAMIC_INIT( core_simple_metrics_per_residue_metrics_CurrentProbabilityMetric )
#endif // SERIALIZATION





