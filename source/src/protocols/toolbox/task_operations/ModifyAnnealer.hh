// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ModifyAnnealer.hh
///
/// @brief Task operation to set high and low temps for annealer as well as whether or not to do a quench step
/// @author Tim Jacobs


#ifndef INCLUDED_protocols_toolbox_task_operations_ModifyAnnealer_HH
#define INCLUDED_protocols_toolbox_task_operations_ModifyAnnealer_HH

// Unit Headers
#include <protocols/toolbox/task_operations/ModifyAnnealer.fwd.hh>

// Core Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>

namespace protocols {
namespace toolbox {
namespace task_operations {

class ModifyAnnealer : public core::pack::task::operation::TaskOperation
{

public:
	ModifyAnnealer();

	ModifyAnnealer(bool disallow_quench, core::Real high_temp, core::Real low_temp);

	virtual ~ModifyAnnealer();

	virtual core::pack::task::operation::TaskOperationOP clone() const;

	virtual void apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const;

	virtual void parse_tag( utility::tag::TagCOP, basic::datacache::DataMap & );
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "ModifyAnnealer"; }

private:

	bool disallow_quench_;
	core::Real high_temp_;
	core::Real low_temp_;

};

} //task_operations
} //toolbox
} //protocols

#endif
