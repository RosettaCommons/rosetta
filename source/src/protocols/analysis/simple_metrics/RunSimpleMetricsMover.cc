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

// Basic/Utility headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <utility>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
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
	protocols::moves::Mover( src ),
	prefix_(src.prefix_),
	suffix_(src.suffix_)
{
	metrics_.clear();
	for ( auto metric: src.metrics_ ) {
		metrics_.push_back( metric ); //Don't need to be copies as we will not be modifying the metrics.
	}
}

RunSimpleMetricsMover::~RunSimpleMetricsMover()= default;

void
RunSimpleMetricsMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	//We should also accept a comma-separated list of previously defined simple metrics
	// Full name of simple_metrics as this is exactly what they are.
	utility::vector1< SimpleMetricCOP > metrics = get_metrics_from_datamap_and_subtags(tag, data);
	set_simple_metrics(metrics);

	prefix_ = tag->getOption< std::string >("prefix", prefix_);
	suffix_ = tag->getOption< std::string >("suffix", suffix_);
}

protocols::moves::MoverOP
RunSimpleMetricsMover::clone() const
{
	return protocols::moves::MoverOP( new RunSimpleMetricsMover( *this ) );
}

protocols::moves::MoverOP
RunSimpleMetricsMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new RunSimpleMetricsMover );
}

// XRW TEMP std::string
// XRW TEMP RunSimpleMetricsMover::get_name() const
// XRW TEMP {
// XRW TEMP  return RunSimpleMetricsMover::mover_name();
// XRW TEMP }

void
RunSimpleMetricsMover::apply( core::pose::Pose & pose, std::string const & prefix, std::string const & suffix){
	set_prefix(prefix);
	set_suffix(suffix);
	apply(pose);
}

void
RunSimpleMetricsMover::apply( core::pose::Pose & pose )
{

	for ( SimpleMetricCOP metric : metrics_ ) {
		debug_assert( metric != nullptr );
		TR << "Running: " << metric->name()<< " - " << "calculating " << metric->metric() << std::endl;
		metric->apply(pose, prefix_, suffix_);

	}
}

void
RunSimpleMetricsMover::add_simple_metric( SimpleMetricCOP metric )
{
	metrics_.push_back( metric );
}

void
RunSimpleMetricsMover::set_simple_metrics( utility::vector1< SimpleMetricCOP > metrics ){
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

void RunSimpleMetricsMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	SimpleMetricFactory::get_instance()->define_simple_metric_xml_schema( xsd );
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "metrics", xs_string, "Comma-separated list of previously defined simple_metrics to be added." )
		+ XMLSchemaAttribute( "prefix", xs_string, "Prefix tag for the values to be added to the output score file for these metrics. (prefix+metric_name+suffix)" )
		+ XMLSchemaAttribute( "suffix", xs_string, "suffix tag for the values to be added to the output score file for these metrics. (prefix+metric_name+suffix)" );
	XMLSchemaSimpleSubelementList subelements;
	subelements.add_group_subelement( & SimpleMetricFactory::get_instance()->simple_metric_xml_schema_group_name );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "Run a set of SimpleMetrics on the pose, with a particular prefix and suffix. (prefix+metric_name+suffix) ", attlist, subelements );
}

std::string RunSimpleMetricsMoverCreator::keyname() const {
	return RunSimpleMetricsMover::mover_name();
}

protocols::moves::MoverOP
RunSimpleMetricsMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new RunSimpleMetricsMover );
}

void RunSimpleMetricsMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RunSimpleMetricsMover::provide_xml_schema( xsd );
}

} // simple_metrics
} // analysis
} // protocols

