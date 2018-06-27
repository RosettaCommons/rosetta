// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/TimingProfileMetric.cc
/// @brief Calculate the time difference between construction and apply/calculate.  Useful to time protocols in RosettaScripts or through mover containers.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <core/simple_metrics/metrics/TimingProfileMetric.hh>
#include <core/simple_metrics/metrics/TimingProfileMetricCreator.hh>

// Core headers
#include <core/simple_metrics/RealMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>

//C++ headers
#include <chrono>

static basic::Tracer TR( "core.simple_metrics.metrics.TimingProfileMetric" );


namespace core {
namespace simple_metrics {
namespace metrics {

using namespace std::chrono;

using minutes_timing = duration< float, std::ratio< 60, 1 >>;
using hours_timing = duration< float, std::ratio< 3600, 1>>;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
TimingProfileMetric::TimingProfileMetric():
	core::simple_metrics::RealMetric()
{
	construction_time_ = high_resolution_clock::now();
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
TimingProfileMetric::~TimingProfileMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
TimingProfileMetric::TimingProfileMetric( TimingProfileMetric const & src ):
	RealMetric( src ),
	construction_time_( src.construction_time_),
	calc_in_hours_( src.calc_in_hours_)
{

}


core::simple_metrics::SimpleMetricOP
TimingProfileMetric::clone() const {
	return SimpleMetricOP(new TimingProfileMetric( *this ) );

}

core::Real
TimingProfileMetric::calculate(const pose::Pose & ) const {
	high_resolution_clock::time_point t2 = high_resolution_clock::now();

	if ( calc_in_hours_ ) {
		return duration_cast< hours_timing >( t2 - construction_time_ ).count();
	} else {
		return duration_cast< minutes_timing >( t2 - construction_time_ ).count();
	}
}

std::string
TimingProfileMetric::name() const {
	return name_static();
}

std::string
TimingProfileMetric::name_static() {
	return "TimingProfileMetric";

}
std::string
TimingProfileMetric::metric() const {
	return "timing_profile";
}

void
TimingProfileMetric::set_calc_in_hours(bool calc_in_hours){
	calc_in_hours_ = calc_in_hours;
}

void
TimingProfileMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &  )
{
	SimpleMetric::parse_base_tag( tag );
	set_calc_in_hours( tag->getOption< bool >( "hours", false));
}

void
TimingProfileMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("hours", xsct_rosetta_bool, "Boolean to set whether we report in hours.  Default is minutes", "false");

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"A metric for measuring a TimingProfile and adding it to the resulting score file.  The time is between construction and calls of apply.  If you use it in RosettaScripts, you can have multiple timings between mover sets and determine the time between them using separate TimingProfileMetrics.   Useful to get runtimes of new movers and applications. ", attlist);
}

void
TimingProfileMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	TimingProfileMetric::provide_xml_schema( xsd );
}

std::string
TimingProfileMetricCreator::keyname() const {
	return TimingProfileMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
TimingProfileMetricCreator::create_simple_metric() const {
	return core::simple_metrics::SimpleMetricOP( new TimingProfileMetric );

}

} //core
} //simple_metrics
} //metrics






