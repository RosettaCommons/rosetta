// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief  Restrict design to residues matching user-specified SASA criteria in the monomeric, bound, or unbound state.
/// @author Jacob Bale (balej@uw.edu)

#ifndef INCLUDED_protocols_toolbox_task_operations_SelectBySASAOperation_hh
#define INCLUDED_protocols_toolbox_task_operations_SelectBySASAOperation_hh

// Unit Headers
#include <protocols/toolbox/task_operations/SelectBySASAOperation.fwd.hh>
#include <protocols/toolbox/task_operations/SelectBySASAOperationCreator.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <string>

// C++ Headers

namespace protocols {
namespace toolbox {
namespace task_operations {

class SelectBySASAOperation : public core::pack::task::operation::TaskOperation {
public:
	SelectBySASAOperation( std::string mode = "sc", std::string state = "monomer", core::Real probe_radius = 2.2, core::Real core_asa = 0, core::Real surface_asa = 30, std::string jump_nums = "1", std::string sym_dof_names = "", bool core = 0, bool boundary = 0, bool surface = 0, bool verbose = 0 );

	virtual ~SelectBySASAOperation();

	virtual core::pack::task::operation::TaskOperationOP clone() const;

	virtual
	void
	apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const;
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "SelectBySASA"; }

	void parse_tag( TagCOP tag , DataMap & );

private:
	std::string mode_;
	std::string state_;
	core::Real probe_radius_;
	core::Real core_asa_;
	core::Real surface_asa_;
	std::string jump_nums_;
	std::string sym_dof_names_;
	bool core_;
	bool boundary_;
	bool surface_;
	bool verbose_;
};

} //namespace task_operations
} //namespace toolbox
} //namespace protocols

#endif
