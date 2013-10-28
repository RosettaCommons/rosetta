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

#include <utility/vector1.hh>


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
	core::Real design_shell() const{ return design_shell_; }
	void include_residue( core::Size const resid );

	virtual ~DesignAroundOperation();

	virtual TaskOperationOP clone() const;

	virtual
	void
	apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const;

	virtual void parse_tag( TagCOP, DataMap & );

	virtual void parse_def( utility::lua::LuaObject const & def );

	void repack_shell( core::Real const repack_shell) {
		repack_shell_ = repack_shell;
		if( repack_shell <= design_shell() )
			design_shell_ = repack_shell;
	}
	core::Real repack_shell() const{ return repack_shell_; }

	void allow_design( bool const a ){ allow_design_ = a; }
	bool allow_design() const{ return allow_design_; }
	void resnums_allow_design( bool const a ){ resnums_allow_design_ = a; }
	bool resnums_allow_design() const{ return resnums_allow_design_; }
	std::set< core::Size > get_designable_residues() const{ return resid_; };
private:
	core::Real design_shell_, repack_shell_; //dflt 8 and 8
	bool allow_design_; //dflt true
	bool resnums_allow_design_; //dflt true
	std::set< core::Size > resid_; // accessible at all times
	std::string string_resnums_; // this can only be parsed at apply time, when the pose is available
};

} //namespace protocols
} //namespace toolbox
} //namespace task_operations

#endif // INCLUDED_protocols_toolbox_TaskOperations_DesignAroundOperation_HH
