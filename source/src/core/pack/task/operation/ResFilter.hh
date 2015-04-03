// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/operation/ResFilter.hh
/// @brief  abstract base class that takes a pose and a residue index, and returns true or false
/// @author ashworth

// do not add any derived classes to this file, unless they are generalized abstract base classes and do not actually 'do any work'

#ifndef INCLUDED_core_pack_task_operation_ResFilter_hh
#define INCLUDED_core_pack_task_operation_ResFilter_hh

// Unit Headers
#include <core/pack/task/operation/ResFilter.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

// C++ Headers

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace task {
namespace operation {

class ResFilter : public utility::pointer::ReferenceCount
{
public:
	typedef pose::Pose Pose;
	typedef utility::tag::TagCOP TagCOP;
public:
	virtual bool operator() ( Pose const &, Size ) const = 0;
	virtual ResFilterOP clone() const = 0;
	/// @brief parser xml-like tags to set class data/parameters
	virtual void parse_tag( TagCOP ) {}
};

// do not add any derived classes to this file, unless they are generalized abstract base classes and do not actually 'do any work'

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core

#endif
