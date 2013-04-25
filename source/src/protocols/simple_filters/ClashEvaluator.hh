// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ClashEvaluator.hh
/// @brief
/// @detailed
///
///
///
/// @author Oliver Lange



#ifndef INCLUDED_protocols_simple_filters_ClashEvaluator_hh
#define INCLUDED_protocols_simple_filters_ClashEvaluator_hh


// Unit Headers

// Package Headers
#include <protocols/evaluation/PoseEvaluator.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <basic/MetricValue.hh>

#include <core/kinematics/Stub.hh>
#include <core/kinematics/RT.hh>

#include <core/io/silent/silent.fwd.hh>


// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

//// C++ headers
#include <list>

namespace protocols {
namespace simple_filters {

template <class T >
class PoseMetricEvaluator : public evaluation::SingleValuePoseEvaluator< T > {
public:
  PoseMetricEvaluator( std::string tag ) : evaluation::SingleValuePoseEvaluator< T >( metric+"_"+tag ) {};
  virtual T apply( core::pose::Pose& pose  ) const;
private:
  std::string metric_;
  std::string tag_;
};


template < class T >
PoseMetricEvaluator::apply(
 pose::Pose& pose
) const {
  basic::MetricValue< T > mr;
  pose.metric(metric_,tag_,mr);
  return mr.value();
}


}
}

#endif
