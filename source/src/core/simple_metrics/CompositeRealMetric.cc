// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/simple_metrics/CompositeRealMetric.cc
///
/// @brief Main class for simple metrics.
/// @author Jared Adolf-Bryfogle ( jadolfbr@gmail.com )

// Unit Headers
#include <core/simple_metrics/CompositeRealMetric.hh>

// Protocol Headers

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/extra_pose_info_util.hh>

namespace core {
namespace simple_metrics {

CompositeRealMetric::CompositeRealMetric():
	SimpleMetric("CompositeRealMetric")
{}


CompositeRealMetric::~CompositeRealMetric() = default;

CompositeRealMetric::CompositeRealMetric( CompositeRealMetric const & src ):
	SimpleMetric( src )
{


}

void
CompositeRealMetric::apply( pose::Pose & pose, std::string prefix, std::string suffix ) const {
	std::map< std::string, core::Real > values = calculate( pose );
	for ( auto value_pair : values ) {
		std::string out_tag = prefix + value_pair.first + "_"+metric() + suffix;
		core::pose::setPoseExtraScore( pose, out_tag, value_pair.second);
	}
}

} //namespace simple_metrics
} //namespace core
