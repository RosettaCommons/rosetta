// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/task_operations/RestrictChainToRepackingOperation.hh
/// @brief  TaskOperation class that restricts a chain to repacking
/// @author Sarel Fleishman sarelf@uw.edu

#ifndef INCLUDED_protocols_toolbox_task_operations_RestrictChainToRepackingOperation_hh
#define INCLUDED_protocols_toolbox_task_operations_RestrictChainToRepackingOperation_hh

// Unit Headers
#include <protocols/toolbox/task_operations/RestrictChainToRepackingOperation.fwd.hh>
#include <protocols/toolbox/task_operations/RestrictOperationsBase.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// Utility Headers
#include <core/types.hh>

// C++ Headers

#include <utility/vector1.hh>


namespace protocols {
namespace toolbox {
namespace task_operations {

/// @details this class is a TaskOperation to prevent repacking of residues not near an interface.
class RestrictChainToRepackingOperation : public RestrictOperationsBase
{
public:
	typedef RestrictOperationsBase parent;

	RestrictChainToRepackingOperation();
	RestrictChainToRepackingOperation( core::Size const chain );
	virtual ~RestrictChainToRepackingOperation();
	virtual TaskOperationOP clone() const;

	virtual
	void
	apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const;

	void chain( core::Size const chain );
	core::Size chain() const;

	virtual void parse_tag( TagCOP, DataMap & );

	static std::string keyname();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	core::Size chain_;
};

} //namespace protocols
} //namespace toolbox
} //namespace task_operations

#endif // INCLUDED_protocols_toolbox_TaskOperations_RestrictChainToRepackingOperation_HH
