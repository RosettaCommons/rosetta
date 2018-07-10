// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/simple_metrics/RealMetric.cc
///
/// @brief Main class for simple metrics.
/// @author Jared Adolf-Bryfogle ( jadolfbr@gmail.com )

// Unit Headers
#include <core/simple_metrics/RealMetric.hh>

// Protocol Headers

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/extra_pose_info_util.hh>

namespace core {
namespace simple_metrics {

RealMetric::RealMetric():
	SimpleMetric("RealMetric")
{}


RealMetric::~RealMetric() = default;

RealMetric::RealMetric( RealMetric const & src ):
	SimpleMetric( src )
{


}

void
RealMetric::apply( pose::Pose & pose, std::string prefix, std::string suffix ) const {

	std::string custom_type = get_custom_type();

	if ( custom_type != "" ) custom_type=custom_type+"_";

	std::string out_tag = prefix + custom_type + metric() + suffix;

	core::Real value = calculate( pose );
	core::pose::setPoseExtraScore( pose, out_tag, value);
}

utility::vector1< std::string >
RealMetric::get_metric_names() const {
	utility::vector1< std::string > names;
	names.push_back( metric() );
	return names;
}

} //namespace simple_metrics
} //namespace core
