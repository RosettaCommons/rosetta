// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/task_operations/SelectResiduesWithinChainOperation.hh
/// @brief  TaskOperation class that restricts a chain to repacking
/// @author Sarel Fleishman sarelf@uw.edu

#ifndef INCLUDED_protocols_toolbox_task_operations_SelectResiduesWithinChainOperation_hh
#define INCLUDED_protocols_toolbox_task_operations_SelectResiduesWithinChainOperation_hh

// Unit Headers
#include <protocols/toolbox/task_operations/SelectResiduesWithinChain.fwd.hh>
#include <protocols/toolbox/task_operations/RestrictOperationsBase.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
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
class SelectResiduesWithinChainOperation : public RestrictOperationsBase
{
public:
	typedef RestrictOperationsBase parent;

	SelectResiduesWithinChainOperation();

	virtual ~SelectResiduesWithinChainOperation();

	virtual TaskOperationOP clone() const;

	virtual
	void
	apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const;
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "SelectResiduesWithinChain"; }

	virtual void parse_tag( TagCOP, DataMap & );

	core::Size chain() const{ return chain_; }
	void chain( core::Size const c ){ chain_ = c; }
	utility::vector1< core::Size > resid() { return resid_; }
	void add_res( core::Size const resid ){ resid_.push_back( resid ); }
	void allow_design( bool const b ){ allow_design_ = b; }
	bool allow_design() const{ return allow_design_; }

	void allow_repacking( bool const b ){ allow_repacking_ = b; }
	bool allow_repacking() const{ return allow_repacking_; }

	void modify_unselected_residues( bool const b ){ modify_unselected_residues_ = b; }
	bool modify_unselected_residues() const { return modify_unselected_residues_; }
private:
	core::Size chain_; // dflt 1;
	utility::vector1< core::Size > resid_; //dflt empty;
	bool allow_design_; //dflt true; allow design at selected positions. If false, restrict to repacking
	bool allow_repacking_; //dflt true; allow repacking in selected residues. If false, prevent repacking in these residues
	bool modify_unselected_residues_; //dflt true; prevent repacking on all residues that are not selected. If false, only affect selected residues; leave unselected residues alone
};

} //namespace protocols
} //namespace toolbox
} //namespace task_operations

#endif // INCLUDED_protocols_toolbox_TaskOperations_SelectResiduesWithinChainOperation_HH
