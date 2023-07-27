// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/analysis/simple_metrics/RunSimpleMetricsMover.cc
/// @brief Run a set of SimpleMetrics (used primarily for RosettaScripts)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/analysis/simple_metrics/RunSimpleMetricsMover.hh>
#include <protocols/analysis/simple_metrics/RunSimpleMetricsMoverCreator.hh>

// Protocol headers
#include <core/simple_metrics/SimpleMetric.hh>
#include <core/simple_metrics/SimpleMetricFactory.hh>
#include <core/simple_metrics/util.hh>


// Core headers
#include <core/pose/Pose.hh>
#include <core/simple_metrics/SimpleMetricData.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>

// Basic/Utility headers
#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.analysis.simple_metrics.RunSimpleMetricsMover" );

namespace protocols {
namespace analysis {
namespace simple_metrics {

using namespace core::simple_metrics;

RunSimpleMetricsMover::RunSimpleMetricsMover():
	protocols::moves::Mover( RunSimpleMetricsMover::mover_name() ),
	metrics_()
{
}

RunSimpleMetricsMover::RunSimpleMetricsMover( utility::vector1< SimpleMetricCOP> const & metrics ):
	protocols::moves::Mover( RunSimpleMetricsMover::mover_name() ),
	metrics_( metrics )
{
}

RunSimpleMetricsMover::RunSimpleMetricsMover( RunSimpleMetricsMover const & src ):
	protocols::moves::Mover( src ) ,
	prefix_(src.prefix_),
	suffix_(src.suffix_),
	metric_to_bfactor_(src.metric_to_bfactor_),
	override_existing_data_(src.override_existing_data_)

{
	metrics_.clear();
	for ( auto const & metric: src.metrics_ ) {
		metrics_.push_back( metric ); //Don't need to be copies as we will not be modifying the metrics.
	}
}

RunSimpleMetricsMover::~RunSimpleMetricsMover()= default;

void
RunSimpleMetricsMover::set_override(bool override_existing_data){
	override_existing_data_ = override_existing_data;
}

/// @brief Provide the citation.
void
RunSimpleMetricsMover::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	using namespace basic::citation_manager;

	citations.add(
		utility::pointer::make_shared< UnpublishedModuleInfo >(
		mover_name(),
		CitedModuleType::Mover,
		"Jared Adolf-Bryfogle",
		"The Scripps Research Institute, La Jolla, CA",
		"jadolfbr@gmail.com"
		)
	);
}

void
RunSimpleMetricsMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) {
	//We should also accept a comma-separated list of previously defined simple metrics
	// Full name of simple_metrics as this is exactly what they are.
	utility::vector1< SimpleMetricCOP > metrics = get_metrics_from_datamap_and_subtags(tag, data);
	set_simple_metrics(metrics);
	set_override(tag->getOption< bool >("override", override_existing_data_));

	prefix_ = tag->getOption< std::string >("prefix", prefix_);
	suffix_ = tag->getOption< std::string >("suffix", suffix_);
	metric_to_bfactor_ = tag->getOption< std::string >("metric_to_bfactor", metric_to_bfactor_);
}

protocols::moves::MoverOP
RunSimpleMetricsMover::clone() const
{
	return utility::pointer::make_shared< RunSimpleMetricsMover >( *this );
}

protocols::moves::MoverOP
RunSimpleMetricsMover::fresh_instance() const
{
	return utility::pointer::make_shared< RunSimpleMetricsMover >();
}


void
RunSimpleMetricsMover::apply( core::pose::Pose & pose, std::string const & prefix, std::string const & suffix, std::string const & metric_to_bfactor){
	set_prefix(prefix);
	set_suffix(suffix);
	set_bfactor(metric_to_bfactor);
	apply(pose);
}

void
RunSimpleMetricsMover::apply( core::pose::Pose & pose )
{
	bool is_custom_type_for_b_factor_valid = false;
	for ( SimpleMetricCOP const & metric : metrics_ ) {
		debug_assert( metric != nullptr );
		TR << "Running: " << metric->name()<< " - " << "calculating " << metric->get_final_sm_type() << std::endl;
		metric->apply(pose, prefix_, suffix_, override_existing_data_);

		// set result of metric to the b factor column if specified by user
		if ( !metric_to_bfactor_.empty() && metric->get_custom_type() == metric_to_bfactor_ ) {
			is_custom_type_for_b_factor_valid = true;
			set_metric_to_bfactor( pose, metric );
		}
	}
	if ( !metric_to_bfactor_.empty() ) {
		runtime_assert_msg( is_custom_type_for_b_factor_valid, "Error in RunSimpleMetrics: Value of metric_to_bfactor does not match any metric custom_type, did you misspell it?");
	}
}

void
RunSimpleMetricsMover::add_simple_metric( SimpleMetricCOP const & metric )
{
	metrics_.push_back( metric );
}

void
RunSimpleMetricsMover::set_simple_metrics( utility::vector1< SimpleMetricCOP > const & metrics ){
	metrics_ = metrics;
}

std::string
RunSimpleMetricsMover::get_name() const {
	return mover_name();
}

std::string
RunSimpleMetricsMover::mover_name() {
	return "RunSimpleMetrics";
}

void
RunSimpleMetricsMover::set_prefix( std::string const & prefix ){
	prefix_ = prefix;
}


void
RunSimpleMetricsMover::set_suffix( std::string const & suffix ){
	suffix_ = suffix;
}

void
RunSimpleMetricsMover::set_bfactor( std::string const & metric_to_bfactor) {
	metric_to_bfactor_ = metric_to_bfactor;
}
void
RunSimpleMetricsMover::set_metric_to_bfactor( core::pose::Pose & pose, SimpleMetricCOP const & metric) {

	core::simple_metrics::SimpleMetricDataOP metric_data = core::simple_metrics::get_sm_data(pose);
	std::map<std::string, std::map<core::Size, core::Real >> all_metric_data_map = metric_data->get_per_residue_real_metric_data();

	std::string const metric_key = prefix_ + metric_to_bfactor_ + "_" + metric->metric() + suffix_;
	auto metric_data_iterator = all_metric_data_map.find(metric_key);

	if ( metric_data_iterator != all_metric_data_map.end() ) {
		std::map<core::Size, core::Real> &metric_data_map = metric_data_iterator->second;
		for ( auto const &entry: metric_data_map ) {
			core::Size const residue_position = entry.first;
			core::Real const residue_metric_value = entry.second;
			for ( core::Size atom_position = 1;
					atom_position <= pose.residue(residue_position).natoms(); ++atom_position ) {
				pose.pdb_info()->temperature(residue_position, atom_position, residue_metric_value);
			}
		}
		TR << "Successfully set the values of " << metric->name() << " as temperature factor." << std::endl;
	} else {
		utility_exit_with_message(
			"Error in RunSimpleMetrics: Value of the metric_to_bfactor does not match any custom_type of a PerResidueRealMetric.");
	}
}



void RunSimpleMetricsMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	SimpleMetricFactory::get_instance()->define_simple_metric_xml_schema( xsd );
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "metrics", xs_string, "Comma-separated list of previously defined simple_metrics to be added." )
		+ XMLSchemaAttribute( "prefix", xs_string, "Prefix tag for the values to be added to the output score file for these metrics. (prefix + custom_type + _ + metric_name + suffix)" )
		+ XMLSchemaAttribute( "suffix", xs_string, "suffix tag for the values to be added to the output score file for these metrics. (prefix + custom_type + _ + metric_name + suffix)" )
		+ XMLSchemaAttribute::attribute_w_default( "override", xsct_rosetta_bool, "Should we override any existing data?", "false")
		+ XMLSchemaAttribute::attribute_w_default( "metric_to_bfactor", xs_string, "Name of a PerResidueRealMetric which values will be written to the b factor column.", "");

	XMLSchemaSimpleSubelementList subelements;
	subelements.add_group_subelement( & SimpleMetricFactory::get_instance()->simple_metric_xml_schema_group_name );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(),
		"Author: Jared Adolf-Bryfogle (jadolfbr@gmail.com)\n"
		"Run a set of SimpleMetrics on the pose, with a particular prefix and suffix. (prefix + custom_type + _ + metric_name + suffix) ", attlist, subelements );
}

std::string RunSimpleMetricsMoverCreator::keyname() const {
	return RunSimpleMetricsMover::mover_name();
}

protocols::moves::MoverOP
RunSimpleMetricsMoverCreator::create_mover() const {
	return utility::pointer::make_shared< RunSimpleMetricsMover >();
}

void RunSimpleMetricsMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RunSimpleMetricsMover::provide_xml_schema( xsd );
}

} // simple_metrics
} // analysis
} // protocols

