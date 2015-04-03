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


#ifndef INCLUDED_protocols_simple_filters_ExternalEvaluator_hh
#define INCLUDED_protocols_simple_filters_ExternalEvaluator_hh


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


//// C++ headers

namespace protocols {
namespace simple_filters {

class ExternalEvaluator : public evaluation::SingleValuePoseEvaluator< core::Real > {
public:
  ExternalEvaluator( std::string tag, std::string command = "echo NO_COMMAND_SPECIFIED" );
	virtual ~ExternalEvaluator();
  virtual core::Real apply( core::pose::Pose& pose ) const;
  virtual bool applicable(  core::pose::Pose const& ) const { return true; };
	void set_command( std::string const&  command_in ) { command_ = command_in; };
	std::string scratch_dir() const {
		return scratch_dir_;
	}
private:
  std::string command_;
	std::string scratch_dir_;
	std::string work_dir_;
};

}
}

#endif
