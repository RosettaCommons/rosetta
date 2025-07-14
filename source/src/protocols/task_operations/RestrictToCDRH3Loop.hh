// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/task_operations/RestrictToCDRH3Loop.hh
/// @brief  This class allows the selection for packing of the Antibody CDR-H3 loop by taking advantage of the PDB numbering schemes that are commonly used for Antibodies
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)

#ifndef INCLUDED_protocols_task_operations_RestrictToCDRH3Loop_hh
#define INCLUDED_protocols_task_operations_RestrictToCDRH3Loop_hh

// Unit headers
#include <protocols/task_operations/RestrictToCDRH3Loop.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// Project headers
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <string>


namespace protocols {
namespace task_operations {

class RestrictToCDRH3Loop : public core::pack::task::operation::TaskOperation
{
public:
	typedef core::pack::task::operation::TaskOperation parent;

public:
	RestrictToCDRH3Loop();

	RestrictToCDRH3Loop(RestrictToCDRH3Loop const & src);

	~RestrictToCDRH3Loop() override;

	core::pack::task::operation::TaskOperationOP clone() const override;

	void
	apply(core::pose::Pose const & pose, core::pack::task::PackerTask & task) const override;
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "RestrictToCDRH3Loop"; }

private:
	bool residue_is_in_h3_loop( core::pose::Pose const & pose, core::Size residue_number ) const;

private:


	// These definitions correspond to the Chothia numbering scheme.  This could be expanded to other numbering schemes in the future
	static core::Size const pdb_numbered_h3_loop_start = 95; // 107;
	static core::Size const pdb_numbered_h3_loop_end = 102;  // 138;
	static std::string const heavy_chain; // = "H";

};

} //namespace task_operations
} //namespace protocols

#endif // include guard
