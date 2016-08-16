// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/task_operations/CrystalContactsOperation.hh
/// @brief  Exclude crystal contacts from design
/// @author Patrick Conway (ptconway@uw.edu)

#ifndef INCLUDED_protocols_toolbox_task_operations_CrystalContactsOperation_hh
#define INCLUDED_protocols_toolbox_task_operations_CrystalContactsOperation_hh

// Unit Headers
#include <protocols/toolbox/task_operations/CrystalContactsOperation.fwd.hh>
#include <protocols/toolbox/task_operations/CrystalContactsOperationCreator.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// Project Headers
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>
#include <core/types.hh>

// Utility Headers
#include <string>

// C++ Headers

namespace protocols {
namespace toolbox {
namespace task_operations {

class CrystalContactsOperation : public core::pack::task::operation::TaskOperation {
public:
	CrystalContactsOperation(
		core::Real all_gap = 0.5,
		core::Real polar_gap = 2.5,
		core::Real max_buried_sasa_ = 0.01,
		bool invert = false,
		bool nbr_radius_to_nbr_radius = false,
		bool nbr_radius_to_atoms = true,
		bool atoms_to_atoms = false);

	virtual ~CrystalContactsOperation();

	virtual core::pack::task::operation::TaskOperationOP clone() const;

	virtual
	void
	apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const;

	void parse_tag( TagCOP tag, DataMap & );
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "CrystalContacts"; }

private:
	bool is_crystal_contact( core::conformation::Residue const & asymm_residue, core::conformation::Residue const & symm_residue ) const;

	core::Real all_gap_;
	core::Real polar_gap_;
	core::Real max_buried_sasa_;
	bool invert_;

	bool nbr_radius_to_nbr_radius_;
	bool nbr_radius_to_atoms_;
	bool atoms_to_atoms_;
};

} //namespace task_operations
} //namespace toolbox
} //namespace protocols

#endif
