// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/composite_metrics/BestMutationsFromProbabilitiesMetric.cc
/// @brief A class for calculating the mutations with the highest delta_probability to the current residues from a PerResidueProbabilitiesMetric.
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)

// Unit headers
#include <core/simple_metrics/composite_metrics/BestMutationsFromProbabilitiesMetric.hh>
#include <core/simple_metrics/composite_metrics/BestMutationsFromProbabilitiesMetricCreator.hh>

// Core headers
#include <core/simple_metrics/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/pointer/memory.hh>

// C++ includes
#include <queue>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "core.simple_metrics.composite_metrics.BestMutationsFromProbabilitiesMetric" );

namespace core {
namespace simple_metrics {
namespace composite_metrics {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
BestMutationsFromProbabilitiesMetric::BestMutationsFromProbabilitiesMetric():
	core::simple_metrics::CompositeRealMetric()
{}

BestMutationsFromProbabilitiesMetric::BestMutationsFromProbabilitiesMetric( core::simple_metrics::PerResidueProbabilitiesMetricCOP metric ):
	core::simple_metrics::CompositeRealMetric()
{
	set_metric( metric );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
BestMutationsFromProbabilitiesMetric::~BestMutationsFromProbabilitiesMetric()= default;

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
BestMutationsFromProbabilitiesMetric::BestMutationsFromProbabilitiesMetric(BestMutationsFromProbabilitiesMetric const &  ) = default;

core::simple_metrics::SimpleMetricOP
BestMutationsFromProbabilitiesMetric::clone() const {
	return utility::pointer::make_shared< BestMutationsFromProbabilitiesMetric >(*this );
}

std::string
BestMutationsFromProbabilitiesMetric::name() const {
	return name_static();
}

std::string
BestMutationsFromProbabilitiesMetric::name_static() {
	return "BestMutationsFromProbabilitiesMetric";

}
std::string
BestMutationsFromProbabilitiesMetric::metric() const {
	return name_static();
}

void
BestMutationsFromProbabilitiesMetric::set_metric(core::simple_metrics::PerResidueProbabilitiesMetricCOP metric){
	metric_ = metric;
}

void
BestMutationsFromProbabilitiesMetric::set_use_cached_data(bool use_cache, std::string const & prefix, std::string const & suffix){
	use_cache_ = use_cache;
	cache_prefix_ = prefix;
	cache_suffix_ = suffix;
	if ( use_cache_ ) {
		TR << "Attempting to use cached data with set prefix/suffix:" <<prefix <<" "<<suffix << std::endl;
	}
}

void
BestMutationsFromProbabilitiesMetric::set_cutoffs( core::Size max_number_mutations, core::Real delta_cutoff ) {
	max_number_mutations_ = max_number_mutations;
	delta_cutoff_ = delta_cutoff;
}

void
BestMutationsFromProbabilitiesMetric::set_fail_on_missing_cache(bool fail){
	fail_on_missing_cache_ = fail;
}

// @brief Compute the ten largest deltas for amino acid mutations.
std::map<std::string, core::Real>
BestMutationsFromProbabilitiesMetric::compute_deltas(
	std::map<core::Size, std::map<core::chemical::AA, core::Real>> const & values,
	core::pose::Pose const & pose
) const
{
	// Priority queue to store largest deltas
	// first a lambda function to define the sorting of the smallest item by probability
	auto comparison = [](std::pair<std::string, core::Real> const & vector_1, std::pair<std::string, core::Real> const & vector_2) {
		return vector_1.second > vector_2.second;
	};
	// then the queue using the function
	std::priority_queue<std::pair<std::string, core::Real>, std::vector<std::pair<std::string, core::Real>>, decltype(comparison)> top_deltas(comparison);

	// loop through all positions
	for ( auto const & position_and_probs : values ) {
		core::Size const & position = position_and_probs.first;
		std::map<core::chemical::AA, core::Real> const & aa_probs = position_and_probs.second;

		// Get the current AA from the pose
		core::chemical::AA current_aa = pose.residue(position).aa();

		// Get its probability
		auto it = aa_probs.find(current_aa);
		core::Real current_prob;
		if ( it != aa_probs.end() ) {
			current_prob = it->second;
		} else {
			current_prob = 0.0;
			TR.Error << "Can't find probability for the current amino acid at position " << position << " setting it to 0.0 instead!" << std::endl;
		}

		// Identify the maximum probability from the other amino acids
		core::Real max_other_prob = 0.0;
		core::chemical::AA max_aa = core::chemical::AA::aa_unk; // to keep compiler happy
		for ( auto const & aa_and_prob : aa_probs ) {
			core::chemical::AA const & aa = aa_and_prob.first;
			core::Real const & prob = aa_and_prob.second;

			if ( aa != current_aa && prob > max_other_prob ) {
				max_other_prob = prob;
				max_aa = aa;
			}
		}

		// Compute the delta and consider for storage
		core::Real delta = max_other_prob - current_prob;
		// only store the delta if it's larger than the cutoff (default=0.0 meaning equally likely as the current amino acid).
		if ( delta >= delta_cutoff_ ) {
			// creating a describing string in "CurrentPositionMutation" format (e.g. "A78V" ) using pose numbering
			std::string mutation_str =
				std::string(1, oneletter_code_from_aa(current_aa)) + std::to_string(position )
				+ std::string(1, oneletter_code_from_aa(max_aa));

			top_deltas.emplace(mutation_str, delta);

			// Maintain number of max mutations
			while ( top_deltas.size() > max_number_mutations_ ) {
				top_deltas.pop();
			}
		}
	}

	// Extract the largest deltas into a map
	std::map<std::string, core::Real> top_deltas_map;
	while ( !top_deltas.empty() ) {
		// save in map
		top_deltas_map.insert( top_deltas.top() );
		top_deltas.pop();
	}

	if ( TR.Info.visible() ) {
		// sort the mutations by position to print them in ascending order
		std::vector<std::pair<int, std::string>> ordered_mutations;
		for ( const auto &kv: top_deltas_map ) {
			const std::string &mutation_str = kv.first;
			std::string position_str = mutation_str.substr(1, mutation_str.size() - 2);
			int position = std::stoi(position_str);
			ordered_mutations.emplace_back(position, mutation_str);
		}

		std::sort(ordered_mutations.begin(), ordered_mutations.end(),
			[](const std::pair<int, std::string> &a, const std::pair<int, std::string> &b) {
				return a.first < b.first;
			}
		);

		TR.Info << "Top mutations, ordered by pose position:" << std::endl;
		for ( const auto &p: ordered_mutations ) {
			const std::string &mutation_str = p.second;
			TR.Info << "Mutation: " << mutation_str << " delta_probability: " << top_deltas_map[mutation_str]
				<< std::endl;
		}
	}
	return top_deltas_map;
}

void
BestMutationsFromProbabilitiesMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data )
{
	SimpleMetric::parse_base_tag( tag );

	core::simple_metrics::SimpleMetricCOP metric = core::simple_metrics::get_metric_from_datamap_and_subtags(tag, data, "metric");

	if ( metric->simple_metric_type() != "PerResidueProbabilitiesMetric" ) {
		utility_exit_with_message("BestMutationsFromProbabilitiesMetric only works with PerResidueProbabilitiesMetrics!");
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

	core::Size max_number_mutations = tag->getOption<core::Size>( "max_mutations", 10);
	auto delta_cutoff = tag->getOption<core::Real>( "delta_cutoff", 0.0 );

	set_cutoffs( max_number_mutations, delta_cutoff );
}

void
BestMutationsFromProbabilitiesMetric::provide_xml_schema(utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;

	attlist + XMLSchemaAttribute::required_attribute( "metric", xs_string, "A PerResidueProbabilitiesMetric to calculate the probability delta between the current amino acids and the most likely ones." );
	attlist + XMLSchemaAttribute::attribute_w_default( "max_mutations", xsct_positive_integer, "The maximum amount of mutations that will be returned (Default=10).", "10" );
	attlist + XMLSchemaAttribute::attribute_w_default( "delta_cutoff", xsct_real, "The cutoff for the delta in probability that will be reported. Default is 0 meaning at least as likely as the current residue", "0.0");
	attlist + XMLSchemaAttribute::attribute_w_default( "use_cached_data",  xsct_rosetta_bool, "Use any data stored in the datacache that matches the set metrics name (and any prefix/suffix.)  Data is stored during a SimpleMetric's apply function, which is called during RunSimpleMetrics", "false");
	attlist + XMLSchemaAttribute("cache_prefix", xs_string, "Any prefix used during apply (RunSimpleMetrics), that we will match on if use_cache is true");
	attlist + XMLSchemaAttribute("cache_suffix", xs_string, "Any suffix used during apply (RunSimpleMetrics), that we will match on if use_cache is true");
	attlist + XMLSchemaAttribute::attribute_w_default("fail_on_missing_cache", xsct_rosetta_bool, "If use_cached_data is True and cache is not found, should we fail?", "true");

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"A CompositeRealMetric for calculating the mutations with the highest delta_probability to the current residues from a PerResidueProbabilitiesMetric. Outputs in the format Mutation-Position-Current, delta_value (e.g. D10A, 0.8)", attlist);

}

std::map< std::string, core::Real >
BestMutationsFromProbabilitiesMetric::calculate(core::pose::Pose const &pose ) const {

	if ( metric_ == nullptr ) {
		utility_exit_with_message("BestMutationsFromProbabilitiesMetric: This metric requires a PerResidueProbabilitiesMetric!");
	}

	std::map< core::Size, std::map< core::chemical::AA, core::Real >> const values = metric_->cached_calculate( pose, use_cache_, cache_prefix_, cache_suffix_, fail_on_missing_cache_ );

	std::map< std::string, core::Real > return_map = compute_deltas( values, pose );

	return return_map;
}

void
BestMutationsFromProbabilitiesMetricCreator::provide_xml_schema(utility::tag::XMLSchemaDefinition & xsd ) const {
	BestMutationsFromProbabilitiesMetric::provide_xml_schema(xsd );
}

std::string
BestMutationsFromProbabilitiesMetricCreator::keyname() const {
	return BestMutationsFromProbabilitiesMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
BestMutationsFromProbabilitiesMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< BestMutationsFromProbabilitiesMetric >();
}

utility::vector1< std::string >
BestMutationsFromProbabilitiesMetric::get_metric_names() const {
	utility::vector1< std::string > names;
	names.push_back( metric() );
	return names;
}

void
BestMutationsFromProbabilitiesMetric::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"BestMutationsFromProbabilitiesMetric", basic::citation_manager::CitedModuleType::SimpleMetric,
		"Moritz Ertelt",
		"University of Leipzig",
		"moritz.ertelt@gmail.com",
		""
		)
	);
}

} //composite_metrics
} //simple_metrics
} //core


#ifdef    SERIALIZATION

template< class Archive >
void
core::simple_metrics::composite_metrics::BestMutationsFromProbabilitiesMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::CompositeRealMetric>( this ) );
	arc( CEREAL_NVP( metric_ ) );
	arc( CEREAL_NVP( use_cache_ ) );
	arc( CEREAL_NVP( cache_prefix_ ) );
	arc( CEREAL_NVP( cache_suffix_ ) );
	arc( CEREAL_NVP( fail_on_missing_cache_ ) );
	arc( CEREAL_NVP ( max_number_mutations_ ) );
	arc( CEREAL_NVP ( delta_cutoff_ ) );

}

template< class Archive >
void
core::simple_metrics::composite_metrics::BestMutationsFromProbabilitiesMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::CompositeRealMetric >( this ) );

	std::shared_ptr< core::simple_metrics::PerResidueProbabilitiesMetric > local_metric;
	arc( local_metric); // CompositeRealMetricCOP
	metric_ = local_metric; // copy the non-const pointer(s) into the const pointer(s)

	arc( metric_ );
	arc( use_cache_ );
	arc( cache_prefix_ );
	arc( cache_suffix_ );
	arc( fail_on_missing_cache_ );
	arc( max_number_mutations_ );
	arc( delta_cutoff_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::simple_metrics::composite_metrics::BestMutationsFromProbabilitiesMetric );
CEREAL_REGISTER_TYPE( core::simple_metrics::composite_metrics::BestMutationsFromProbabilitiesMetric )

CEREAL_REGISTER_DYNAMIC_INIT( core_simple_metrics_composite_metrics_BestMutationsFromProbabilitiesMetric )
#endif // SERIALIZATION





