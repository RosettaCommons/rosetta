// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/simple_metrics/CompositeStringMetric.cc
///
/// @brief Main class for simple metrics.
/// @author Jared Adolf-Bryfogle ( jadolfbr@gmail.com )

// Unit Headers
#include <core/simple_metrics/CompositeStringMetric.hh>

// Protocol Headers

// Core headers
#include <core/pose/Pose.hh>
#include <core/simple_metrics/util.hh>
#include <core/simple_metrics/SimpleMetricData.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace simple_metrics {

CompositeStringMetric::CompositeStringMetric():
	SimpleMetric("CompositeStringMetric")
{}


CompositeStringMetric::~CompositeStringMetric() = default;

CompositeStringMetric::CompositeStringMetric( CompositeStringMetric const & src ):
	SimpleMetric( src )
{


}

std::map< std::string, std::string>
CompositeStringMetric::cached_calculate(pose::Pose const & pose, bool use_cache, std::string prefix, std::string suffix, bool fail_on_missing_cache) const {
	std::string name = prefix + get_final_sm_type() + suffix;

	if ( use_cache && has_sm_data( pose ) ) {
		std::map< std::string, std::string > value;
		bool data_found = get_sm_data(pose)->get_value(name, value);
		if ( data_found ) {
			return value;
		} else if ( fail_on_missing_cache ) {
			utility_exit_with_message("Could not find CompositeStringMetric: "+name+" in pose");
		} else {
			return calculate(pose);
		}
	} else {
		return calculate(pose);
	}
}

void
CompositeStringMetric::apply( pose::Pose & pose, std::string prefix, std::string suffix ) const {

	std::string out_tag = prefix + get_final_sm_type() + suffix;
	std::map< std::string, std::string > values = calculate( pose );
	MetricKey mk;
	get_sm_data(pose)->set_value( mk, out_tag, values );

	//for ( auto value_pair : values ) {
	// std::string out_tag = prefix + value_pair.first + "_" + get_final_sm_type() + suffix;
	// core::pose::setPoseExtraScore( pose, out_tag, value_pair.second);
	//}
}

} //namespace simple_metrics
} //namespace core

#ifdef    SERIALIZATION



template< class Archive >
void
core::simple_metrics::CompositeStringMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::SimpleMetric>( this ) );

}

template< class Archive >
void
core::simple_metrics::CompositeStringMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::SimpleMetric >( this ) );

}

SAVE_AND_LOAD_SERIALIZABLE( core::simple_metrics::CompositeStringMetric );
CEREAL_REGISTER_TYPE( core::simple_metrics::CompositeStringMetric )

CEREAL_REGISTER_DYNAMIC_INIT( core_simple_metrics_CompositeStringMetric )
#endif // SERIALIZATION

