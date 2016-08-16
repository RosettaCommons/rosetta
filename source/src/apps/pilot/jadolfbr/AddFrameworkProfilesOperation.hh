// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /AddFrameworkProfilesOperation.hh
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_task_operations_AddFrameworkProfilesOperation_hh
#define INCLUDED_protocols_antibody_task_operations_AddFrameworkProfilesOperation_hh

#include <protocols/antibody/task_operations/AddFrameworkProfilesOperation.fwd.hh>
#include <protocols/antibody/AntibodyInfo.fwd.hh>
#include <protocols/antibody/AntibodyEnum.hh>

#include <core/pack/task/operation/TaskOperation.hh>

// Utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace antibody {
namespace task_operations {


/// @brief Add Framework Profiles as the task operation for a framework chain.
/// See protocols/toolbox/task_operations/ResidueProbTaskOperation for more.
///
class AddFrameworkProfilesOperation : public core::pack::task::operation::TaskOperation {
public:

	AddFrameworkProfilesOperation();
	AddFrameworkProfilesOperation(AntibodyInfoCOP ab_info);
	AddFrameworkProfilesOperation(AntibodyInfoCOP ab_info, bool use_conservative_set);

	virtual ~AddFrameworkProfilesOperation();

	AddFrameworkProfilesOperation(AddFrameworkProfilesOperation const & src);

	virtual
	void
	apply(core::pose::Pose const & pose, core::pack::task::PackerTask & task) const;

	/// @brief Use the set of conservative mutations instead of profiles.
	void
	set_use_conservative_set( bool use_conservatives);

private:
	AntibodyInfoCOP ab_info_;
	//std::string species_;
	//std::string gene_;
	//std::string chain_;

	//Profile Options
	core::Size picking_rounds_;
	bool add_to_current_;
	bool use_native_type_;
	bool use_conservative_set_;

	///Needed for default and RS constructor.
	AntibodyNumberingSchemeEnum numbering_scheme_;
	CDRDefinitionEnum cdr_definition_;

};

} //task_operations
} //antibody
} //protocols

#endif	//INCLUDED_protocols_antibody_task_operations_AddFrameworkProfilesOperation_hh



