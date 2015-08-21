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


#ifndef INCLUDED_protocols_sparta_ChemicalShiftEvaluator_hh
#define INCLUDED_protocols_sparta_ChemicalShiftEvaluator_hh


// Unit Headers

// Package Headers
//#include <protocols/evaluation/ExternalEvaluator.hh>
#include <protocols/evaluation/PoseEvaluator.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

#include <protocols/sparta/Sparta.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


//// C++ headers

namespace protocols {
namespace sparta {

class ChemicalShiftEvaluator : public protocols::evaluation::SingleValuePoseEvaluator< core::Real > {
public:
	ChemicalShiftEvaluator( std::string tag, std::string cst_file );
	virtual core::Real apply( core::pose::Pose& pose ) const;
	virtual bool applicable( core::pose::Pose const&pose ) const;
private:
	mutable sparta::Sparta sparta_; //since this class is far from const-correct..
};

// class ChemicalShiftEvaluator : public ExternalEvaluator {
// public:
//   ChemicalShiftEvaluator( std::string tag, std::string cst_file );
//   //  virtual core::Real apply( core::pose::Pose& pose ) const;
//   virtual bool applicable( core::pose::Pose const&pose ) const;
// private:

// };

}
}

#endif
