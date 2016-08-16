// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ClashEvaluator.hh
/// @brief
/// @details
///
///
///
/// @author Oliver Lange


#ifndef INCLUDED_protocols_simple_filters_PoseMetricEvaluator_hh
#define INCLUDED_protocols_simple_filters_PoseMetricEvaluator_hh


// Unit Headers

// Package Headers

// Project Headers
#include <core/pose/Pose.hh> // replace this with .fwd by moving methods to .cc file!
#include <basic/MetricValue.hh>


// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers

#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace simple_filters {

template <class T >
class PoseMetricEvaluator : public evaluation::SingleValuePoseEvaluator< T > {
public:
	PoseMetricEvaluator( std::string const& metric, std::string const& tag ) :
		evaluation::SingleValuePoseEvaluator< T >( metric+"_"+tag ),
		metric_( metric ),
		tag_( tag )
	{};
	virtual T apply( core::pose::Pose& pose  ) const;
private:
	std::string metric_;
	std::string tag_;
};


template < class T >
T PoseMetricEvaluator<T >::apply(
	core::pose::Pose& pose
) const {
	basic::MetricValue< T > mr;
	pose.metric(metric_,tag_,mr);
	return mr.value();
}


}
}

#endif
