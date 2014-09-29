// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @detailed
///	  Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange



#ifndef INCLUDED_protocols_simple_filters_MinValueEvaluator_hh
#define INCLUDED_protocols_simple_filters_MinValueEvaluator_hh


// Unit Headers
#include <protocols/evaluation/PoseEvaluator.hh>

// Package Headers

// Project Headers
#include <core/io/silent/silent.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

//// C++ headers
#include <list>

// namespace protocols {
// namespace simple_filters {

// template< class T >
// class MinValueEvaluator : public evaluation::SingleValuePoseEvaluator<T> {
// public:
//   MinValueEvaluator( utility::pointer::owning_ptr< evaluation::SingleValuePoseEvaluator<T> const > eval_in ) :
// 		evaluation::SingleValuePoseEvaluator< T > ( "min_"+eval_in->name() ),
// 		evaluator_( eval_in ),
// 		min_value_( 1000000 )
//   {};

//   virtual
//   T apply( core::pose::Pose& pose ) const {
//     T val = evaluator_->apply( pose );
//     if ( min_value_ > val ) min_value_ = val;
//     return min_value_;
//   }

// private:
//   utility::pointer::owning_ptr< evaluation::SingleValuePoseEvaluator<T> const > evaluator_;
//   mutable T min_value_;
// };

// }
// }

#endif
