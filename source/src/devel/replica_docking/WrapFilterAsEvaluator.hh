// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @details
/// @author Zhe Zhang


#ifndef INCLUDED_devel_replica_docking_WrapFilterAsEvaluator_hh
#define INCLUDED_devel_replica_docking_WrapFilterAsEvaluator_hh


// Unit Headers

// Package Headers
//#include <protocols/evaluation/ExternalEvaluator.hh>
#include <protocols/evaluation/PoseEvaluator.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

#include <protocols/filters/Filter.hh>
#include <protocols/filters/Filter.fwd.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


//// C++ headers

namespace devel {
namespace replica_docking {

class WrapFilterAsEvaluator : public protocols::evaluation::SingleValuePoseEvaluator< core::Real > {
public:
	WrapFilterAsEvaluator( protocols::filters::FilterOP, std::string );
	core::Real apply( core::pose::Pose& pose ) const override;
	bool applicable( core::pose::Pose const&pose ) const override;
private:
	protocols::filters::FilterOP filter_; //since this class is far from const-correct..
};

// class WrapFilterAsEvaluator : public ExternalEvaluator {
// public:
//   WrapFilterAsEvaluator( std::string tag, std::string cst_file );
//   //  virtual core::Real apply( core::pose::Pose& pose ) const;
//   virtual bool applicable( core::pose::Pose const&pose ) const;
// private:

// };

}
}

#endif
