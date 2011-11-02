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



#ifndef INCLUDED_protocols_evaluation_RPF_ScoreEvaluator_hh
#define INCLUDED_protocols_evaluation_RPF_ScoreEvaluator_hh

//#include <protocols/evaluation/RPF_ScoreEvaluator.fwd.hh>


// Unit Headers
#include <protocols/evaluation/PoseEvaluator.hh>

// Package Headers
// AUTO-REMOVED #include <protocols/noesy_assign/NoesyModule.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>


// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <protocols/noesy_assign/CrossPeakList.fwd.hh>
#include <utility/vector1.hh>


// AUTO-REMOVED #include <utility/vector1.hh>

//// C++ headers
// AUTO-REMOVED #include <list>

namespace protocols {
namespace evaluation {

class RPF_ScoreEvaluator : public evaluation::SingleValuePoseEvaluator< core::Real > {
public:
  RPF_ScoreEvaluator( std::string tag, core::Real dcut = 5 );

	virtual core::Real apply( core::pose::Pose& pose ) const;
	virtual bool applicable(  core::pose::Pose const& pose ) const;

private:
	mutable noesy_assign::CrossPeakListOP crosspeaks_;
	core::Real dcut_;
};

}
}

#endif
