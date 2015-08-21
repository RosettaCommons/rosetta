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
/// @details
///
///
///
/// @author Oliver Lange


#ifndef INCLUDED_protocols_constraints_additional_CombinedConstraintEvaluator_hh
#define INCLUDED_protocols_constraints_additional_CombinedConstraintEvaluator_hh

#ifdef WIN32
#include <protocols/constraints_additional/ConstraintEvaluator.hh>
#endif


// Unit Headers

// Package Headers
#include <protocols/evaluation/PoseEvaluator.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

#ifdef WIN32 // for visual studio
#else
#endif

#include <core/io/silent/silent.fwd.hh>


// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <protocols/constraints_additional/ConstraintEvaluator.fwd.hh>
#include <utility/vector1.hh>


//// C++ headers

namespace protocols {
namespace constraints_additional {

class CombinedConstraintEvaluator : public evaluation::PoseEvaluator {
public:
	CombinedConstraintEvaluator( std::string tag, std::string filename, Size constraints_combine_ratio_ = 2, Size repeat = 10 );
	//sets xxx_cst and xxx_viol columns
	virtual void apply( core::pose::Pose& pose, std::string tag, core::io::silent::SilentStruct &pss) const;

	using evaluation::PoseEvaluator::apply;

	//returns constraint score
	virtual core::Real apply( core::pose::Pose& pose ) const;

	virtual core::Size size() const { return 2; }
	virtual std::string name( core::Size i ) const;

private:
	utility::vector1< ConstraintEvaluator > cst_lib_;
	std::string name_;

};

}
}

#endif
