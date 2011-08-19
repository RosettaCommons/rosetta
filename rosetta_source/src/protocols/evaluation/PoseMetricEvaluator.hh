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



#ifndef INCLUDED_protocols_evaluation_PoseMetricEvaluator_hh
#define INCLUDED_protocols_evaluation_PoseMetricEvaluator_hh


// Unit Headers

// Package Headers
// AUTO-REMOVED #include <protocols/evaluation/PoseEvaluator.hh>

// Project Headers
#include <core/pose/Pose.hh> // replace this with .fwd by moving methods to .cc file!
#include <basic/MetricValue.hh>

// AUTO-REMOVED #include <core/kinematics/Stub.hh>
// AUTO-REMOVED #include <core/kinematics/RT.hh>

// AUTO-REMOVED #include <core/io/silent/silent.fwd.hh>


// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
// AUTO-REMOVED #include <utility/vector1.hh>

//// C++ headers
// AUTO-REMOVED #include <list>

//Auto Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>


namespace protocols {
namespace evaluation {

template <class T >
class PoseMetricEvaluator : public evaluation::SingleValuePoseEvaluator< T > {
public:
  PoseMetricEvaluator( std::string const& metric, std::string const& tag ) :
    SingleValuePoseEvaluator< T >( metric+"_"+tag ),
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
