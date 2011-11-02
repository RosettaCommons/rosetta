// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file JumpEvaluator.hh
/// @brief
/// @detailed
///
///
///
/// @author Oliver Lange



#ifndef INCLUDED_protocols_evaluation_BatchNrEvaluator_hh
#define INCLUDED_protocols_evaluation_BatchNrEvaluator_hh


// Unit Headers

// Package Headers
#include <protocols/evaluation/PoseEvaluator.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

#include <core/io/silent/silent.fwd.hh>


// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


// AUTO-REMOVED #include <utility/vector1.hh>

//// C++ headers
// AUTO-REMOVED #include <list>

namespace protocols {
namespace evaluation {

//@brief yields a column with the number of the batch decoy was evaluated in
class BatchEvaluator : public evaluation::PoseEvaluator {
public:
  virtual void apply( core::pose::Pose& pose, std::string tag, core::io::silent::SilentStruct &pss ) const;

	virtual core::Size size() const { return 1; };
	virtual std::string name( core::Size ) const { return "batch"; };
private:
};



//@brief yields a column with the number of the batch decoy was evaluated in
class BatchNrEvaluator : public evaluation::SingleValuePoseEvaluator< core::Size > {
public:
  BatchNrEvaluator() : SingleValuePoseEvaluator< core::Size >( "batch_nr" ) {};
  virtual core::Size apply( core::pose::Pose& pose  ) const;
private:
};


}
}

#endif
