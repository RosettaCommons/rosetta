// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/task_operations/RestrictAAsFromProbabilities.cc
/// @brief A class to restrict designable amino acids depending on the probabilities from a PerResidueProbabilitiesMetric
/// @author Moritz Ertelt (moritz.ertelt@gmail.com)

#include <protocols/task_operations/RestrictAAsFromProbabilities.hh>
#include <protocols/task_operations/RestrictAAsFromProbabilitiesCreator.hh>

// core headers
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/task_op_schemas.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/simple_metrics/util.hh>
#include <core/pack/task/ResidueLevelTask.hh>
#include <core/chemical/ResidueTypeSet.hh>

// basic
#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/pointer/memory.hh>

static basic::Tracer TR( "protocols.task_operations.RestrictAAsFromProbabilities" );

namespace protocols {
namespace task_operations {
using namespace core::pack::task::operation;

RestrictAAsFromProbabilities::RestrictAAsFromProbabilities():
	TaskOperation()

{

}

RestrictAAsFromProbabilities::~RestrictAAsFromProbabilities() {}

TaskOperationOP
RestrictAAsFromProbabilities::clone() const {
	return TaskOperationOP( new RestrictAAsFromProbabilities( *this ) );
}

void
RestrictAAsFromProbabilities::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap& data ){

	// PerResidueProbabilitiesMetric we will sample from
	core::simple_metrics::SimpleMetricCOP metric = core::simple_metrics::get_metric_from_datamap_and_subtags(tag, data, "metric");

	if ( metric->simple_metric_type() != "PerResidueProbabilitiesMetric" ) {
		utility_exit_with_message("The RestrictAAsFromProbabilities TaskOperation only works with PerResidueProbabilitiesMetrics! Your metric " + metric->name() + " is of type " + metric->simple_metric_type() );
	}

	core::simple_metrics::PerResidueProbabilitiesMetricCOP res_metric = utility::pointer::dynamic_pointer_cast< core::simple_metrics::PerResidueProbabilitiesMetric const>( metric );

	set_metric( res_metric );

	set_prob_cutoff( tag->getOption< core::Real >( "prob_cutoff", prob_cutoff_ ) );
	set_delta_prob_cutoff( tag->getOption< core::Real>( "delta_prob_cutoff", delta_prob_cutoff_ ) );

	// options for using cached data
	bool use_cache = tag->getOption< bool >("use_cached_data", false);
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
RestrictAAsFromProbabilities::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task) const {

	using namespace core::chemical;

	// check that it's a PerResidueProbabilitiesMetric
	if ( metric_ == nullptr ) {
		utility_exit_with_message("RestrictAAsFromProbabilities: This TaskOperation requires a PerResidueProbabilitiesMetric!");
	}

	// get values from the metric
	TR << "Calculating/Fetching probabilities from metric..." << std::endl;
	std::map< core::Size, std::map< AA, core::Real >> values = metric_->cached_calculate( pose, use_cache_, cache_prefix_, cache_suffix_, fail_on_missing_cache_ );
	TR << " Done!" << std::endl;
	// disable AAs that don't match the user set probability cutoffs
	disable_AAs( values, pose, task );
	TR << "Disabled AAs with overall probability lower than " << prob_cutoff_ << " and delta to current lower than  " << delta_prob_cutoff_ << " ." << std::endl;
}

void
RestrictAAsFromProbabilities::set_use_cached_data(bool use_cache, std::string const & prefix, std::string const & suffix){
	use_cache_ = use_cache;
	cache_prefix_ = prefix;
	cache_suffix_ = suffix;
	if ( use_cache_ ) {
		TR << "Attempting to use cached data with set prefix/suffix:" <<prefix <<" "<<suffix << std::endl;
	}
}
void
RestrictAAsFromProbabilities::set_fail_on_missing_cache(bool fail){
	fail_on_missing_cache_ = fail;
}
void
RestrictAAsFromProbabilities::set_metric(core::simple_metrics::PerResidueProbabilitiesMetricCOP metric){
	metric_ = metric;
}
void RestrictAAsFromProbabilities::set_prob_cutoff(core::Real prob_cutoff) {
	prob_cutoff_ = prob_cutoff;
}
void RestrictAAsFromProbabilities::set_delta_prob_cutoff(core::Real delta_prob_cutoff) {
	delta_prob_cutoff_ = delta_prob_cutoff;
}

void
RestrictAAsFromProbabilities::disable_AAs(
	std::map<core::Size, std::map<core::chemical::AA, core::Real>>& probabilities,
	core::pose::Pose const& pose,
	core::pack::task::PackerTask& task
) const {
	using namespace core::chemical;

	for ( auto & pos_and_probs : probabilities ) {
		core::Size position = pos_and_probs.first;
		std::map<core::chemical::AA, core::Real>& aa_probs = pos_and_probs.second;

		AA current_aa = pose.residue( position ).aa();

		core::Real current_prob;
		if ( aa_probs.find(current_aa) != aa_probs.end() ) {
			current_prob = aa_probs.at( current_aa );
		} else {
			current_prob = 0.0;
			TR.Error << "The amino acid type at position " << position << " was not found in the probabilities map!" << std::endl;
		}

		core::chemical::ResidueTypeSetCOP residue_set = pose.residue_type_set_for_pose();
		utility::vector1< std::string > allowed_base_names;

		// Loop through the amino acid probabilities and set the allowed amino acids
		for ( auto aa_probabilitiy : aa_probs ) {
			AA aa = aa_probabilitiy.first;
			core::Real prob = aa_probabilitiy.second;

			if ( prob >= prob_cutoff_ && ( prob - current_prob ) >= delta_prob_cutoff_ ) {
				core::chemical::ResidueTypeCOP res_type = residue_set->get_representative_type_aa( aa );
				allowed_base_names.push_back( res_type->base_name() );
				TR.Debug << "At position " << position << " we allowed " << res_type->base_name() << " with prob " << prob << std::endl;
			}
		}
		task.nonconst_residue_task( position ).restrict_restypes( allowed_base_names );
	}
}


std::string
RestrictAAsFromProbabilities::keyname() {
	return "RestrictAAsFromProbabilities";
}

void
RestrictAAsFromProbabilities::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute( "metric", xs_string, "A PerResidueProbabilitiesMetric to restrict amino acid identities from." );
	attlist + XMLSchemaAttribute::attribute_w_default(
		"prob_cutoff", xsct_real,
		"Probability cutoff for amino acids to consider for sampling. Defaults to 0.0001", "0.0001");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"delta_prob_cutoff", xsct_real,
		"Delta probability (current-proposed) cutoff for amino acids to consider for sampling. Defaults to 0 (as likely as the current)", "0");

	attlist + XMLSchemaAttribute::attribute_w_default( "use_cached_data",  xsct_rosetta_bool, "Use any data stored in the datacache that matches the set metrics name (and any prefix/suffix.)  Data is stored during a SimpleMetric's apply function, which is called during RunSimpleMetrics", "false");
	attlist + XMLSchemaAttribute("cache_prefix", xs_string, "Any prefix used during apply (RunSimpleMetrics), that we will match on if use_cache is true");
	attlist + XMLSchemaAttribute("cache_suffix", xs_string, "Any suffix used during apply (RunSimpleMetrics), that we will match on if use_cache is true");
	attlist + XMLSchemaAttribute::attribute_w_default("fail_on_missing_cache", xsct_rosetta_bool, "If use_cached_data is True and cache is not found, should we fail?", "true");

	task_op_schema_w_attributes( xsd, keyname(), attlist, "A TaskOperation to restrict designable amino acids depending on the probabilities from a PerResidueProbabilitiesMetric" );
}

TaskOperationOP
RestrictAAsFromProbabilitiesCreator::create_task_operation() const
{
	return utility::pointer::make_shared< RestrictAAsFromProbabilities >();
}

std::string
RestrictAAsFromProbabilitiesCreator::keyname() const
{
	return RestrictAAsFromProbabilities::keyname();
}

void
RestrictAAsFromProbabilitiesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RestrictAAsFromProbabilities::provide_xml_schema( xsd );
}


} //task_operations
} //protocols
