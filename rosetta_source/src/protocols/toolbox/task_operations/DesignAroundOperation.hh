// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/DesignAroundOperation.hh
/// @brief  TaskOperation class that restricts a chain to repacking
/// @author Sarel Fleishman sarelf@uw.edu

#ifndef INCLUDED_protocols_toolbox_task_operations_DesignAroundOperation_hh
#define INCLUDED_protocols_toolbox_task_operations_DesignAroundOperation_hh

// Unit Headers
#include <protocols/toolbox/task_operations/DesignAroundOperation.fwd.hh>
#include <protocols/toolbox/task_operations/RestrictOperationsBase.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

// Utility Headers
#include <core/types.hh>

// C++ Headers
#include <string>
#include <set>

namespace protocols {
namespace toolbox {
namespace task_operations {

///@details this class is a TaskOperation to prevent repacking of residues not near an interface.
class DesignAroundOperation : public RestrictOperationsBase
{
public:
	typedef RestrictOperationsBase parent;

	DesignAroundOperation();

	void design_shell( core::Real const radius );
	void include_residue( core::Size const resid );

	virtual ~DesignAroundOperation();

	virtual TaskOperationOP clone() const;

	virtual
	void
	apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const;

	virtual void parse_tag( TagPtr );

	bool repack_on() const;
	void repack_on( bool const repack_on );

private:
	core::Real design_shell_;
	std::set< core::Size > resid_; // accessible at all times
	std::string string_resnums_; // this can only be parsed at apply time, when the pose is available
	bool repack_on_; // dflt true; do we leave all non-designed positions as repack or prevent repacking?
};

} //namespace protocols
} //namespace toolbox
} //namespace task_operations

#endif // INCLUDED_protocols_toolbox_TaskOperations_DesignAroundOperation_HH
